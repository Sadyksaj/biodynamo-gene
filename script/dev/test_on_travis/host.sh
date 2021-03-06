#!/bin/bash
# Runs a build on a docker container emulating travis-ci.
# Arguments:
#  $1 Branch name that should be built

if [[ $# -ne 1 ]]; then
  echo "Wrong number of arguments.
Usage:
Runs a build on a docker container emulating travis-ci.
Arguments:
  $1 Branch name that should be built"
  exit 1
fi

BRANCH=$1

# enables GUI apps
xhost +local:root

sudo docker stop travis-14
sudo docker rm travis-14

sudo docker run --name travis-14 --net=host --env="DISPLAY" -v //var/run/docker.sock:/var/run/docker.sock -dit travisci/ci-garnet:packer-1496954857 /sbin/init

# get path of this script
pushd `dirname $0` > /dev/null
SCRIPTPATH=`pwd`
popd > /dev/null

sudo docker cp $SCRIPTPATH/inside_docker.sh travis-14:/
sudo docker exec -it travis-14 chmod +x /inside_docker.sh
sudo docker exec -it travis-14 ./inside_docker.sh $BRANCH

RETURN_VAL=$?
if [ "$RETURN_VAL" == "0" ]; then
  echo "Test successful"
  exit $RETURN_VAL
else
  echo "Test FAILED"
fi
