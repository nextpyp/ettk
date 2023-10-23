# Directory containing class headers.
SET(VTK_GEOVIS_HEADER_DIR "/dscrhome/ab690/code/ETTK/dcc/external/VTK5.10.1/Geovis")

# Classes in vtkGeovis.
SET(VTK_GEOVIS_CLASSES
  "vtkCompassRepresentation"
  "vtkCompassWidget"
  "vtkGeoAdaptiveArcs"
  "vtkGeoAlignedImageRepresentation"
  "vtkGeoAlignedImageSource"
  "vtkGeoArcs"
  "vtkGeoAssignCoordinates"
  "vtkGeoCamera"
  "vtkGeoFileImageSource"
  "vtkGeoFileTerrainSource"
  "vtkGeoGlobeSource"
  "vtkGeoGraticule"
  "vtkGeoImageNode"
  "vtkGeoInteractorStyle"
  "vtkGeoProjectionSource"
  "vtkGeoProjection"
  "vtkGeoRandomGraphSource"
  "vtkGeoSampleArcs"
  "vtkGeoSource"
  "vtkGeoSphereTransform"
  "vtkGeoTerrain"
  "vtkGeoTerrain2D"
  "vtkGeoTerrainNode"
  "vtkGeoTransform"
  "vtkGeoTreeNode"
  "vtkGeoTreeNodeCache"
  "vtkGeoView"
  "vtkGeoView2D"
  "vtkGlobeSource")

# Abstract classes in vtkGeovis.
SET(VTK_GEOVIS_CLASSES_ABSTRACT
  "vtkGeoSource")

# Wrap-exclude classes in vtkGeovis.
SET(VTK_GEOVIS_CLASSES_WRAP_EXCLUDE)

# Wrap-special classes in vtkGeovis.
SET(VTK_GEOVIS_CLASSES_WRAP_SPECIAL)

# Wrappable non-class headers for vtkGeovis.
SET(VTK_GEOVIS_WRAP_HEADERS)

# Set convenient variables to test each class.
FOREACH(class ${VTK_GEOVIS_CLASSES})
  SET(VTK_CLASS_EXISTS_${class} 1)
ENDFOREACH(class)
FOREACH(class ${VTK_GEOVIS_CLASSES_ABSTRACT})
  SET(VTK_CLASS_ABSTRACT_${class} 1)
ENDFOREACH(class)
FOREACH(class ${VTK_GEOVIS_CLASSES_WRAP_EXCLUDE})
  SET(VTK_CLASS_WRAP_EXCLUDE_${class} 1)
ENDFOREACH(class)
FOREACH(class ${VTK_GEOVIS_CLASSES_WRAP_SPECIAL})
  SET(VTK_CLASS_WRAP_SPECIAL_${class} 1)
ENDFOREACH(class)
