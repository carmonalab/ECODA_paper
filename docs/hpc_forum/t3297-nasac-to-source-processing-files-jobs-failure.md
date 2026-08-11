# NASAC to source processing files:jobs failure

- Source: https://hpc-community.unige.ch/t/3297

- Created: 2024-02-08T16:59:16.123Z

- Posts: 2

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Nada.Kojovic (2024-02-08T16:59:16.227Z)

## Primary informations
Username: KOJOVIC
Cluster: BAOBAB
## Description
I tried to preprocess a few videos stored on Nasac. The job starts as usual but fails without any log being recorded.
To mount nasac i used:
dbus-launch bash
/usr/bin/gio mount smb://nasac-m2.unige.ch/m-AutismLab
To launch the processing:
python run_openpose.py --input_videos /home/kojovic/.gvfs/smb-share:server=nasac-m2.unige.ch,share=m-autismlab/AutismVideos/OP/OPTEST
Below I enclose the content of the 2 files needed for this step, namely run_openpose.py and template_2020.sh.
Thank you so much for your help!
Nada
```
##### run_openpose.py#####
#!/usr/bin/env python
from __future__ import print_function
import sys
import argparse
import os
import glob
import re
import uuid
import shutil
from subprocess import call

parser = argparse.ArgumentParser(description='Run openpose on a folder of images')
parser.add_argument('--input_videos', metavar='I', type=str, required=True,
                    help='Path to folder containing videos')
### ATTENTION, these models could not correspond to the OpenPose
### version present in the Docker image, better to stick with the
### latter
##### ### ATTENTION, they are not automatically loaded from the canonical
##### ### location, thus copying them :'(
##### parser.add_argument('--openpose_models_dir', type=str, default="/home/kojovic/openpose-master/models", help="Path to Openpose models dir")

args = parser.parse_args()
input_videos = os.path.abspath(args.input_videos)
##### openpose_models_dir = os.path.abspath(args.openpose_models_dir)

print("="*80)
print("Reading Videos...")
videos = glob.glob(os.path.join(input_videos, "*.mp4"))
print(videos)
print("Found {0} mp4 videos inside {1}".format(len(videos), input_videos))

print("Validating videos...")
for _idx, _video in enumerate(videos):
	_video = _video.split("/")
	source_filename = _video[-1].strip()
	target_filename = _video[-1].strip()
	if "_".join(target_filename.split()) != target_filename:
		print("Found whitespace in filename {0}, replacing it with '_'".format(target_filename))
		target_filename = "_".join(target_filename.split())
		source_path = "/".join(_video[:-1]+[source_filename])
		target_path = "/".join(_video[:-1]+[target_filename])
		os.rename(source_path, target_path)
		videos[_idx] = target_path

print("Videos validate !")
print("Creating JOBS folder...")
JOBS_PATH = os.path.join(input_videos, "JOBS")
try:
	shutil.rmtree(JOBS_PATH)
except:
	pass
try:
	os.makedirs(JOBS_PATH)
except:
	pass
print("JOBS FOLDER created at : {0}".format(JOBS_PATH))

print("Creating OUTPUT Folder....")
OUTPUT_PATH = os.path.join(input_videos, "OUTPUT")
try:
	shutil.rmtree(OUTPUT_PATH)
except:
	pass
try:
	os.makedirs(OUTPUT_PATH)
except:
	pass
print("OUTPUT FOLDER created at : {0}".format(OUTPUT_PATH))

print("Creating CONFIG folder...")
CONFIG_PATH = os.path.join(input_videos, "CONFIG")
try:
        shutil.rmtree(CONFIG_PATH)
except:
        pass
try:
        os.makedirs(CONFIG_PATH)
except:
        pass
print("CONFIG FOLDER created at : {0}".format(CONFIG_PATH))
print("Creating LOGS folder...")
LOGS_PATH = os.path.join(input_videos, "LOGS")
try:
        shutil.rmtree(LOGS_PATH)
except:
        pass
try:
        os.makedirs(LOGS_PATH)
except:
        pass
print("LOGS FOLDER created at : {0}".format(LOGS_PATH))

for _idx, _video in enumerate(videos):
	input_filename  = _video.split("/")[-1]
	template_2020_BLACK = open("template_2020_BLACK.sh").read()
	template_2020_BLACK = template_2020_BLACK.replace("INSERT_LOGS_DIRECTORY_HERE", LOGS_PATH)
	template_2020_BLACK = template_2020_BLACK.replace("INSERT_INPUT_DIRECTORY_HERE", input_videos)
	template_2020_BLACK = template_2020_BLACK.replace("INSERT_OUTPUT_DIRECTORY_HERE", OUTPUT_PATH)
	template_2020_BLACK = template_2020_BLACK.replace("INSERT_FILENAME_HERE", input_filename)
	target_filepath = os.path.join(CONFIG_PATH, str(uuid.uuid4())+".sh")
	print("Writing config file for {0} at {1}".format(input_filename, target_filepath))
	fp = open(target_filepath, "w")
	fp.write(template_2020_BLACK)
	fp.close()

os.chdir(LOGS_PATH)
### ATTENTION, see above
##### print("Copying over Open Pose models directory")
##### shutil.copytree(openpose_models_dir, "./models")
# Run jobs
for _config in glob.glob(os.path.join(CONFIG_PATH, "*.sh")):
	print("PUSHING JOB : ", _config)
	call(["sbatch", _config])
```
```
#### template_2020.sh ####

#!/bin/bash

#SBATCH --output=%J-test_adf.log
#SBATCH --error=%J-test_adf.err
#SBATCH --partition=shared-gpu
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --gres=gpu:2

### <https://hpc-community.unige.ch/t/tutorial-launch-openpose-with-gpu-support-through-singularity/593>
module load GCCcore/8.2.0 Singularity/3.4.0-Go-1.12
echo "Starting Processing at `date`   "
cd INSERT_LOGS_DIRECTORY_HERE

FILENAME="INSERT_FILENAME_HERE"
INPUT_DIR="INSERT_INPUT_DIRECTORY_HERE"
OUTPUT_DIR="INSERT_OUTPUT_DIRECTORY_HERE"
OUTPUT_FOLDER=$OUTPUT_DIR/$FILENAME

echo $INPUT_DIR
echo $FILENAME

echo "Creating output folder ${OUTPUT_FOLDER}..."
rm -rf $OUTPUT_FOLDER
mkdir $OUTPUT_FOLDER
mkdir $OUTPUT_FOLDER/json

### ATTENTION, the original Docker image is available at
### <https://hub.docker.com/r/garyfeng/docker-openpose/>
###   NB, the "--pwd /openpose-master" option is mandatory since the
###       build/examples/openpose/openpose.bin is compiled with
###       relative-path dymanic libraries...
DOCKER_OPENPOSE_IMG_PATH="${HOME}/openpose_orchestrator/openpose.simg"
DOCKER_OPENPOSE_IMG_WORKING_DIRECTORY='/openpose-master'
DOCKER_OPENPOSE_IMG_OPENPOSE_BIN='build/examples/openpose/openpose.bin'
CMD="\
srun singularity exec --nv --pwd ${DOCKER_OPENPOSE_IMG_WORKING_DIRECTORY} ${DOCKER_OPENPOSE_IMG_PATH} ${DOCKER_OPENPOSE_IMG_OPENPOSE_BIN} \
--no_display \
--net_resolution '-1x368' \
--scale_number 4 \
--scale_gap 0.25 \
--video $INPUT_DIR/$FILENAME \
--write_coco_json $OUTPUT_FOLDER/json/keypoints_coco.json \
--write_video $OUTPUT_FOLDER/result.avi"

echo "Running the following command:"
echo "$CMD"
echo
time \
 eval ${CMD}

echo "JOB Complete at `date`   "

# Convert result.avi to result.mp4
#ffmpeg -i $OUTPUT_FOLDER/result.avi  -strict 2 $OUTPUT_FOLDER/result.mp4
#rm -rf $OUTPUT_FOLDER/result.avi

exit
```


## Post 2 by @Yann.Sagon (2024-03-04T09:42:35.129Z)

Hi  @Nada.Kojovic[@Nada.Kojovic](https://hpc-community.unige.ch/u/nada.kojovic)
First of all, I’m not sure it is possible to access the NASAC files from the singularity container if you mount it outside of the container.
You said you mounted the NASAC using:
```
dbus-launch bash
/usr/bin/gio mount smb://nasac-m2.unige.ch/m-AutismLab
```
Are you able to read the files (and write if necessary)?
The files are only available on the node where you issued the command. It means if you did it on the login node, the files won’t be available on a compute node where you job will run. You should add those lines in your sbatch to access the NASAC files from the compute nodes.
Regarding my first comment: if it isn’t working, you probably need to mount the NASAC inside the container.
A remark about your sbatch: do not use module load to use Singularity, it is available directly since a while.
Best
Yann
