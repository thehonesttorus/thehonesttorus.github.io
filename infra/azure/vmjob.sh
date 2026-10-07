#!/bin/bash
# Runs on the experiment VM (invoked through an Azure managed run command; WORKSAS arrives as a protected parameter).
# Environment: JOB=<name>, CMD=<shell command, run inside /opt/work/src>, WORKSAS=<container SAS URL>.
# Downloads work/bundle/codebundle.tgz into /opt/work/src, runs CMD in the background with OUT=out/$JOB,
# then uploads out/$JOB (log.txt, DONE with the exit code, and anything CMD wrote there) to work/results/$JOB/.
set -u
: "${JOB:?}" "${CMD:?}" "${WORKSAS:?}"
BASE="${WORKSAS%%\?*}"; TOKEN="${WORKSAS#*\?}"
for i in $(seq 120); do [ -f /opt/work/READY ] && break; sleep 10; done   # cloud-init installs venv and azcopy first
[ -f /opt/work/READY ] || { echo "cloud-init not finished after 20 min"; exit 1; }
mkdir -p /opt/work/src && cd /opt/work/src
azcopy copy "$BASE/bundle/codebundle.tgz?$TOKEN" /opt/work/codebundle.tgz --log-level ERROR >/dev/null && tar xzf /opt/work/codebundle.tgz -C /opt/work/src
mkdir -p "out/$JOB"
printf '%s\n' "$CMD" > "/opt/work/cmd_$JOB.sh"
cat > "/opt/work/run_$JOB.sh" <<INNER
#!/bin/bash
cd /opt/work/src
export PATH=/opt/venv/bin:\$PATH OUT=out/$JOB
bash "/opt/work/cmd_$JOB.sh" > "out/$JOB/log.txt" 2>&1
echo \$? > "out/$JOB/DONE"
azcopy copy "out/$JOB/*" "$BASE/results/$JOB/?$TOKEN" --recursive --log-level ERROR >/dev/null
INNER
chmod +x "/opt/work/run_$JOB.sh"
nohup "/opt/work/run_$JOB.sh" >/dev/null 2>&1 &
echo "launched $JOB (pid $!) on $(hostname), $(nproc) cores"
