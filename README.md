# Soft-Auto-HMC+RWMH tests scripts

## Usage

### Clusters

The `nextflow.config` file instructs Nextflow to use Apptainer with a custom Docker image, so the user only needs to run e.g.
```bash
./nextflow-sockeye.sh run funnel_scale.nf
```
for Sockeye or
```bash
./nextflow-alliance.sh run funnel_scale.nf
```
for the [Alliance](https://alliancecan.ca/en) clusters.


### Local

If you have Julia and cmdstan installed locally, you can run
```bash
./nextflow run funnel_scale.nf
```
**Note**: this assumes that there is an environment variable `CMDSTAN` pointing to the directory where cmdstan lives. You can do this via
```bash
export CMDSTAN=/full/path/to/cmdstan-vx.xx
```

## Recommended `.bashrc` settings for Sockeye

Add the following at the bottom of your `.bashrc`, replacing

- `YOUR_ALLOC_NAME` with the name of the allocation you are assigned to, and
- `YOUR_GITHUB_KEY` with the filename of your private ssh key for accessing Github.

```bash
# export the name of the allocation you belong to
export ALLOC_NAME="YOUR_ALLOC_NAME"

# load modules
module load gcc/9.4.0
module load openjdk
module load git
module load apptainer

# start ssh-agent if not already started, then add keys
if [ -z "$SSH_AUTH_SOCK" ] ; then
  eval `ssh-agent -s` > /dev/null
  ssh-add $HOME/.ssh/YOUR_GITHUB_KEY >/dev/null 2>&1
fi
```
