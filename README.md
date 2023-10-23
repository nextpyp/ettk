# Electron Tomography ToolKit (ETTK)

Main components:
- `virus_segment_membrane`: segmentation of closed surfaces using min-cut/max-flow routines
- `virus_segment_membrane_select_threshold`: visualization of implicit surfaces to assess segmentation quality
- `Correlation3DNew`: geometric template search constrained to a surface
- `CutVolumes3DFromPositions`: extraction of sub-volumes from tomograms
- `LoopCreateVolumeList`: creation of metadata for sub-volume averaging
- `Test_Metric_Filter`: filtering and masking of sub-volumes
- `MPI_Classification`: sub-volume alignment, averaging and classification using spherical harmonics

# Compilation instructions

## 1. Go to the desired working directory and save that location into the MYPATH variable
```
export MYPATH=`pwd`
```

## 2. Clone the ETTK repository
```
git clone git@github.com:nextpyp/ettk.git
cd ettk
```

## 3. Build ettk-devel.sif container if needed (sudo access required)
`sudo singularity build ettk-devel.sif singularity/ettk-devel.def`

## 4. Enter the development environment (bind paths with -B if needed, e.g., -B /work)
`singularity shell ettk-devel.sif`

## 5. Build VTK
```
mkdir ${MYPATH}/ettk/external/VTK5.10.1/build
cd ${MYPATH}/ettk/external/VTK5.10.1/build
ccmake -Wno-dev -D CMAKE_CXX_COMPILER=/usr/bin/g++ -D CMAKE_C_COMPILER=/usr/bin/gcc ..
	Type "c" for configure
	Type "c" for configure
	Type "g" to generate
make -j 4
```

## 6. Build ITK
```
mkdir ${MYPATH}/ettk/external/InsightToolkit-4.2.1/build
cd ${MYPATH}/ettk/external/InsightToolkit-4.2.1/build
ccmake -Wno-dev -D CMAKE_CXX_COMPILER=/usr/bin/g++ -D CMAKE_C_COMPILER=/usr/bin/gcc ..
	Type "c" for configure
	Type "c" for configure
	Type "g" to generate
make -j 4
```

## 7. Build binaries
```
cd ${MYPATH}/ettk
rm CMakeCache.txt*
ccmake -D ITK_DIR=${MYPATH}/ettk/external/InsightToolkit-4.2.1/build -D VTK_DIR=${MYPATH}/ettk/external/VTK5.10.1/build -D BLITZ_DIR=${MYPATH}/ettk/external/blitz-0.9 -D CMAKE_CXX_COMPILER=/usr/bin/g++ -D CMAKE_C_COMPILER=/usr/bin/gcc .
	Type "c" for configure
	Type "c" for configure
	Type "g" to generate
make -j 4
```
