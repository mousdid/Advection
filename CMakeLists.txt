##################
# Début de l'en-tête.
# Normalement, vous n'aurez jamais rien à modifier dans cette partie.
cmake_minimum_required( VERSION 3.0 )

set( CXX_STANDARD_REQUIRED On )
set( CXX_STANDARD C++11 )

set( LIBRARY_PATH_EIGEN "../../../libraries/eigen/" )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "build/" )

# Spécifie les options passées au compilateur. Ici : l'optimisation du code binaire et les warnings.
set( CMAKE_CXX_FLAGS "-Wall" )
# set( CMAKE_CXX_FLAGS_DEBUG "-g" ) # Les paramètres par défaut: décommentez si vous voulez les changer
# set( CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" ) # Les paramètres par défaut: décommentez si vous voulez les changer

# Options diverses de Cmake, à ne pas modifier.
set( CMAKE_SKIP_INSTALL_RULES True )

# Fin de l'en-tête
##################

##################
# Options spécifiques à ce projet.

project( PG201_TP5 VERSION 2023 )

# Décommentez la ligne que vous souhaitez, puis recompilez.
#set( CMAKE_BUILD_TYPE "Debug" )
 set( CMAKE_BUILD_TYPE "Release" )

# Nom de l'exécutable qui sera créé dans le répertoire "build/".
# Nous ajoutons "Debug" ou "Release" au nom de l'exécutable afin de facilement les différencier.
set( TARGET_NAME "tp5_fv_2d_${CMAKE_BUILD_TYPE}" )

# Liste *exhaustive* des fichiers à compiler ensemble.
set( PROJECT_COMPILATION_FILE_LIST src/libraries/Data/DataFile.cpp src/libraries/Data/Function.cpp src/libraries/Mesh/Mesh2D.cpp src/FiniteVolume.cpp src/TimeScheme.cpp src/main.cpp )

# Précise que l'exécutable sera à assembler avec ces fichiers compilés.
add_executable( ${TARGET_NAME} ${PROJECT_COMPILATION_FILE_LIST} )

# Ajoute les répertoires des bibliothèques liées, ici Eigen et les petites biblio locales au projet.
target_include_directories( ${TARGET_NAME} PUBLIC "${LIBRARY_PATH_EIGEN}" "src/libraries" )

# Fin des options spécifiques à ce projet.
##################

#################
## Rappels et notes
#
# - Pour générer le fichier 'Makefile', lancez "cmake CMakeLists.txt".
# - Pour compiler votre code, lancez simplement "make".
# - Pour voir en détail les commandes de compilations qu'exécute "make", lancez "make VERBOSE=1".
# - Pour nettoyer votre répertoire de compilation, lancez "make clean".
#################
