#!/bin/bash
# Package a bundle exactly as for upload, validate the archive, extract it, and re-run the
# EXTRACTED estimator (client/server emulation, first 2 dev MLPs + the 1024x32 smoke shape via the
# subprocess runner) so the numbers can be compared with the unpackaged runs.
# usage: verify_package.sh BUNDLE_NAME   (bundles/<name>/ -> packages/<name>.tar.gz)
set -e
N=$1
HERE=$(cd "$(dirname "$0")/.." && pwd)
export PATH=/root/whest/bin:$PATH
mkdir -p $HERE/packages $HERE/results/package
OUT=$HERE/packages/$N.tar.gz
rm -f $OUT
(cd $HERE/bundles/$N && rm -rf __pycache__ && whest package --estimator . --output $OUT --yes > $HERE/results/package/${N}_package.log 2>&1)
whest validate-package $OUT > $HERE/results/package/${N}_validate.log 2>&1 && echo "$N: validate-package OK"
X=$(mktemp -d /tmp/claude-0/sub/pkgx-XXXX)
tar xzf $OUT -C $X
tar tzvf $OUT > $HERE/results/package/${N}_contents.txt
( cd $X && sha256sum estimator.py ) | sed "s#  #  packaged:#" > $HERE/results/package/${N}_sha.txt
( cd $HERE/bundles/$N && sha256sum estimator.py ) | sed "s#  #  bundle:#" >> $HERE/results/package/${N}_sha.txt
cat $HERE/results/package/${N}_sha.txt
/root/whest/bin/python $HERE/scripts/grader_emul.py run --estimator $X/estimator.py --dataset /tmp/claude-0/sub/dev2 \
    --tag ${N}_packaged_emul --out $HERE/results/package/${N}_packaged_emul_dev2.json --url ipc:///tmp/fs-pkg.sock
/root/whest/bin/python $HERE/scripts/measure_run.py --estimator $X/estimator.py --dataset /tmp/claude-0/sub/rob_w1024_d32 \
    --runner subprocess --max-threads 2 --tag ${N}_packaged_d32 --out $HERE/results/package/${N}_packaged_sub_d32.json
