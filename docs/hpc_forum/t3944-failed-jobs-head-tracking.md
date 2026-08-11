# Failed jobs head tracking

- Source: https://hpc-community.unige.ch/t/3944

- Created: 2025-05-01T22:16:15.242Z

- Posts: 1

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Nada.Kojovic (2025-05-01T22:16:15.285Z)

## Primary informations
Username: $kojovic
Cluster: $BAOBAB
## Description
If i launch the below provided script it gets cancelled.
## Steps to Reproduce
run sbatch tracking_extraction.job
actual script:
```
#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --partition=shared-gpu
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --array=0-9  # If you have 10 videos
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

## Activate virtual environment
source /home/share/schaer2/anshul/miniconda3/etc/profile.d/conda.sh
conda activate video_extraction

# Pick the correct line from video_list.txt based on SLURM_ARRAY_TASK_ID
# Each sub-job processes one video.
VIDEO=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" /srv/beegfs/scratch/shares/schaer_sensit/ESCS/video_list.txt)
OUT_PATH="/home/share/schaer2/anshul/tracks"

echo "Processing video: $VIDEO"

# Call your processing script with the selected VIDEO
python 1_track_people.py --video-file $VIDEO --output-dir $OUT_PATH
```
## Expected Result
iterate over videos stored in : /srv/beegfs/scratch/shares/schaer_sensit/ESCS/video_list.txt
## Actual Result
image
image1080×172 82.3 KB
[image1080×172 82.3 KB](https://hpc-community.unige.ch//hpc-community.unige.ch/uploads/default/original/1X/5b040aca550ec99d5f12695ee964fe4068f0a6e3.jpeg)
