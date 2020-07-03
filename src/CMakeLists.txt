#Generate a ROOT dictionary for vector<vector<unisgned int>>

# Geerate custom dictionaries
ROOT_GENERATE_DICTIONARY(  CustomClassesDict
                            ../inc/Units.h
                            ../inc/Utils.h
                            ../inc/Run.h
                            ../inc/ Pmt.h
                           ../inc/Waveform.h
                           ../inc/ChargeFit.h
                           ../inc/HighChargeFit.h
	                         LINKDEF
                           ../inc/CustomClasses_LinkDef.h
                         )

if (${ROOT_VERSION} VERSION_GREATER "6.0")
	install(
		FILES
		${CMAKE_CURRENT_BINARY_DIR}/libCustomClassesDict_rdict.pcm
		${CMAKE_CURRENT_BINARY_DIR}/libCustomClassesDict.rootmap
		DESTINATION lib/)
endif (${ROOT_VERSION} VERSION_GREATER "6.0")

#Create the CustomClasses library with the ROOT dictionary.
add_library(  CustomClasses
              SHARED
              Utils.cc
              Run.cc
              Pmt.cc
              Waveform.cc
              ChargeFit.cc
              HighChargeFit.cc
            )

#Link to ROOT.
target_link_libraries(CustomClasses ${ROOT_LIBRARIES})
