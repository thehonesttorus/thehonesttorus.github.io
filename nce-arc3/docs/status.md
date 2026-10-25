# Status: blind scoreboard and count-only diagnostics

Governing rule: no frame was looked at, no game source read. Everything below is either the
harness scoreboard (levels completed, action counts, RHAE) or an agent-internal count.

## How the diagnostics were read

`tools/harness.py` prints, per game, the levels completed, the action counts per level and the
agent's `diagnostics()` dictionary: how many demonstrations per action id, how many exact
hypotheses per action (`H_sizes`), how many behaviour classes synthesis could not fit
(`unexplained_classes`), how many combinations failed and by how many cells
(`synth_mismatch_hist`), how many clock colours the level has, step times, and how often each
decision rule fired (try-once, planning, Kato, fallthrough). Census scripts in the scratchpad
aggregated the same kind of counts over stored transitions (changed-cell histograms, how many
components moved/reshaped/appeared per record, whether the pointed cell changed, whether a move
sometimes failed with identical local surroundings).

## What the counts said, and what was changed because of them

| Observation (counts only) | Theory reading | Change |
|---|---|---|
| Nearly every transition in three games changed one to four cells inside a larger same-colour component, never the pointed cell's component as a whole | The pointed selector was coarser than the atom the action acts on; the cell is the finest projection | `SelPointed` gained its cell member; sub-component recolour classes |
| A small two-colour region changed on every action of every kind | The action-independent flow σ_t, not a correspondence of any action | Clock colours: a colour whose post-state is a function of its own pre-state, witnessed under two different actions or data; the exact test is taken in the quotient by its cells |
| One action was chosen 56 times in a row while being a no-op | Partial words made an action perpetually "newborn"; ties went to the first candidate | Newborn is relative to each word's own context; certified no-ops and deaths are never re-chosen; ties go to the least-demonstrated action |
| A component's shape changed in place with a congruent D4 image overlapping its old cells | Two readings in Parry order: identity with a changed rendering (cost 0) before a Rot letter (cost log 7) | Readings in synthesis |
| A move sometimes did not happen with identical colours at the destination | The gate is not local to the mover's surroundings; a history or hidden state | Collision gates (component, colour, inside) were added; the remaining cases are open |
| Cached perceptions drifted (colour 0 became 15) | The engine mutates its grid in place; the perception cache aliased it | Frames are copied on construction |
| Synthesis exhausted its budget on games with dozens of components | Gate fitting re-evaluated the same predicate on the same component thousands of times | Per-record gate tables (lazy), superset pre-filter, positives-first predicate filter, a transport cache per algebra |
| A step took seconds on frames with two hundred components | The matcher built its cost matrix pair by pair in Python; animation sub-steps were perceived and decomposed though no rule reads them | Vectorised cost matrix; perception memoised per frame and rule set; animation sub-steps kept as raw frame pairs. Step time fell by roughly a factor of three |
| Whether one action's move happened was separated exactly by the mass of one colour | A gate on an occupation number of the moment map (a frame-level observable) | Meter gates: a colour mass against a threshold read from the demonstrations, for move clauses and death gates (letter `Meter`, dimension 2·16·4097) |
| A colour shrank by one or two cells per click, the count depending on the click | Powers of a boundary letter, gated by the datum | `Erode^k`, `Grow^k`; gates that test the pointed component (`ptr`) |
| Two colours were taken for a clock although a null click left them unchanged | A clock is a change no action withholds, with a witness of action independence | Null steps count against a clock; the same own-frame change must be seen under two action contexts |

## Scoreboard

(filled by the latest `tools/harness.py` run over the 25 public games; see the table appended
below by the run script)
