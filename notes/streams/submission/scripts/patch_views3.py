"""Residual cut, round 3 (applied on top of patch_views2.py): per-key invalidation of the
view memo.

    python patch_views3.py bundles/v29r2/estimator.py bundles/v29r3/estimator.py

Rounds 1-2 cleared the whole memo whenever ANY Strassen pool entry was allocated or grown
(133 times in a worker's first predict, 37 in its second), so the first two MLPs of every
worker re-sliced almost as much as upstream V29. Now:
- the first allocation of a new pool key invalidates nothing (nothing can reference it);
- growing key K drops only the memo entries that view K's old buffer (tracked by owner key),
  so the old buffer is released and every other cached view survives.
Same ops, same order, same data: FLOPs and outputs must be bit-identical.
"""
import sys

REPL = []
REPL.append(('''        self._sc = {}   # residual cut: cached slices of pooled buffers
        self._vm = {}   # residual cut: id(view) -> (view, {tag: derived views}), pooled views only
        self._ext = {}  # residual cut (round 2): id(permanent pool buffer) -> (buffer, {key: view})

    def _reset_views(self):
        self._sc.clear()
        self._vm.clear()

    def _reg(self, o):
        self._vm[id(o)] = (o, {})
        return o

    def _derived(self, X, tag, make):
        e = self._vm.get(id(X))
        if e is None or e[0] is not X:
            return make()
        r = e[1].get(tag)
        if r is None:
            r = make()
            if isinstance(r, tuple):
                for o in r:
                    self._reg(o)
            else:
                self._reg(r)
            e[1][tag] = r
        return r
''', '''        self._sc = {}   # residual cut: cached slices of pooled buffers
        self._vm = {}   # residual cut: id(view) -> (view, {tag: derived views}, owner pool key)
        self._ext = {}  # residual cut (round 2): id(permanent pool buffer) -> (buffer, {key: view})
        self._owned = {}  # round 3: pool key -> ids of memo entries viewing that key's buffer
        self._sck = {}    # round 3: pool key -> its _sc keys

    def _purge(self, key):
        """Round 3: pool `key` was reallocated -- forget only the views of its old buffer."""
        for i in self._owned.pop(key, ()):
            e = self._vm.get(i)
            if e is not None and e[2] == key:
                del self._vm[i]
        for k in self._sck.pop(key, ()):
            self._sc.pop(k, None)

    def _reg(self, o, owner=None):
        self._vm[id(o)] = (o, {}, owner)
        if owner is not None:
            self._owned.setdefault(owner, []).append(id(o))
        return o

    def _derived(self, X, tag, make):
        e = self._vm.get(id(X))
        if e is None or e[0] is not X:
            return make()
        r = e[1].get(tag)
        if r is None:
            r = make()
            if isinstance(r, tuple):
                for o in r:
                    self._reg(o, e[2])
            else:
                self._reg(r, e[2])
            e[1][tag] = r
        return r
'''))
REPL.append(('''        k = (key, sel)
        s = self._sc.get(k)
        if s is None or s[0] is not buf:
            s = (buf, self._reg(buf[0] if sel == 0 else buf[:sel]))
            self._sc[k] = s
        return s[1]
''', '''        k = (key, sel)
        s = self._sc.get(k)
        if s is None or s[0] is not buf:
            s = (buf, self._reg(buf[0] if sel == 0 else buf[:sel], key))
            self._sc[k] = s
            self._sck.setdefault(key, []).append(k)
        return s[1]
'''))
REPL.append(('''        b = self.pool.get(key)
        if b is None or b.shape[0] < shape[0]:
            b = fnp.empty(shape, dtype=self.dtype)
            fnp.copyto(b, 0.0)   # metered page touch (first-MLP residual, F80)
            self.pool[key] = b
            self._reset_views()
        return b
''', '''        b = self.pool.get(key)
        if b is None or b.shape[0] < shape[0]:
            if b is not None:
                self._purge(key)
            b = fnp.empty(shape, dtype=self.dtype)
            fnp.copyto(b, 0.0)   # metered page touch (first-MLP residual, F80)
            self.pool[key] = b
        return b
'''))
REPL.append(('''        if e is None or e[0].shape[0] < b:
            b5 = fnp.empty(shape5, dtype=self.dtype)
            fnp.copyto(b5, 0.0)   # metered page touch (first-MLP residual, F80)
            b4 = fnp.reshape(b5, (b, shape5[1] * shape5[2], shape5[3], shape5[4]))
            e = (b5, b4)
            self.pool[key] = e
            self._reset_views()
        k = ("b2", key, b)
        s = self._sc.get(k)
        if s is None or s[0] is not e[0]:
            s = (e[0], self._reg(e[0][:b]), self._reg(e[1][:b]))
            self._sc[k] = s
        return s[1], s[2]
''', '''        if e is None or e[0].shape[0] < b:
            if e is not None:
                self._purge(key)
            b5 = fnp.empty(shape5, dtype=self.dtype)
            fnp.copyto(b5, 0.0)   # metered page touch (first-MLP residual, F80)
            b4 = fnp.reshape(b5, (b, shape5[1] * shape5[2], shape5[3], shape5[4]))
            e = (b5, b4)
            self.pool[key] = e
        k = ("b2", key, b)
        s = self._sc.get(k)
        if s is None or s[0] is not e[0]:
            s = (e[0], self._reg(e[0][:b], key), self._reg(e[1][:b], key))
            self._sc[k] = s
            self._sck.setdefault(key, []).append(k)
        return s[1], s[2]
'''))


def main():
    src, dst = sys.argv[1:3]
    t = open(src).read()
    for old, new in REPL:
        assert t.count(old) == 1, old[:160]
        t = t.replace(old, new)
    assert "_reset_views" not in t
    import os
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    open(dst, "w").write(t)


if __name__ == "__main__":
    main()
