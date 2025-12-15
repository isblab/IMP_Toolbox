# Scripts to check sampling quality

### housekeeping
Run the script as follows:

```bash
# for single trajectory
python housekeeping.py -i /path/to/trajectory_1 -r single -b 0.9
# for multiple trajectories
python housekeeping.py -i /path/to/all_trajectories/ -r multiple -b 0.9
```

- `-i` : path to trajectory directory (single) or it's super-directory (multiple)
- `-r` : run_type (single or multiple)
- `-b` : burn_in fraction of frames to discard

Additional arguments:
- `-c` : num_proc i.e. number of processes to use (default: 4)
- `--restraints`: list of restraints to plot the score trajectories for (e.g. ["ExcludedVolumeSphere:EVR", "GaussianEMRestraint:EMR"])
- `--detect_equilibration`: whether to detect score equilibration

### MCMC acceptance rates
Modify `test_acceptance_ball_mover.py` .

### Seeing swap rates, score equilibration for different restraints, getting auto-correlation time
Modify `housekeeping.py`
