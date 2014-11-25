#!/bin/bash

# VERSION CHECK
ROOT_DIR=..
VERSION_FILE="$ROOT_DIR/Version.txt"
if [ ! -f $VERSION_FILE ]; then
    echo "File not found! : $VERSION_FILE"
    exit 1
fi

# TODO: CLEAN BEFORE PACKAGING
##############################

# TEMPORARY FOLDER CREATION
VERSION_DATE=`date +"%Y-%m-%d-%H-%M-%S"`
GUDHI="_GUDHI_"
VERSION_REVISION=`cat $VERSION_FILE`
VERSION_DIR="$VERSION_DATE$GUDHI$VERSION_REVISION"
echo $VERSION_DIR
mkdir "$VERSION_DIR"

# TOP LEVEL FILE COPY
cp $VERSION_FILE $VERSION_DIR
cp $ROOT_DIR/README $VERSION_DIR
cp $ROOT_DIR/Conventions.txt $VERSION_DIR
cp $ROOT_DIR/COPYING $VERSION_DIR

# PACKAGE LEVEL COPY
PACKAGE_INC_DIR="/include"
PACKAGE_SRC_DIR="/source"
PACKAGE_EX_DIR="/example"
for package in $ROOT_DIR/src/*
do
  echo $package
  if [ -d "$package" ]
  then
    if [ "$package" == "common" ]
    then
      cp -r $package $VERSION_DIR
    fi
    if [ -d "$package$PACKAGE_INC_DIR" ]
    then
      cp -r $package$PACKAGE_INC_DIR $VERSION_DIR
    fi
    if [ -d "$package$PACKAGE_SRC_DIR" ]
    then
      cp -r $package$PACKAGE_SRC_DIR $VERSION_DIR
    fi
    if [ -d "$package$PACKAGE_EX_DIR" ]
    then
      cp -r $package$PACKAGE_EX_DIR $VERSION_DIR
    fi
  fi
done

# ZIP DIR AND REMOVE IT
tar -zcf "$VERSION_DIR.tar.gz" "$VERSION_DIR"
rm -rf "$VERSION_DIR"
