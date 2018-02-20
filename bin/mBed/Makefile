#
# File: Makefile
#
# Form:	make
# Path:	D572/Gordon/Clustal-Adapt
# Note:	adaptation of Clustal code to Gordons SeedMap
# Meta:	Clustal, SeedMap
# Date:	2009-03-26, Thur
#	2009-04-01, Wed, moved from MakeGordon to Makefile
#

CC = gcc
CXX = g++
FAB_FLAGS= #-DDIAGNOSTIC 

COPT     = -O2
CFLAGS   = -g $(FAB_FLAGS) $(COPT) -DOS_UNIX
CXXFLAGS = -g $(FAB_FLAGS) $(COPT) -DOS_UNIX
LIBDIR   = -L. 
LIBS     = -lClustal
LIBSS	 = -lSparse -lClustal
INCLUDES = -I../ -I./ -I./pairwise -I./substitutionMatrix -I./multipleAlign -I./general -I./interface -I./fileInput -I./alignment -I./tree -I./tree/UPGMA -I./sparseMap

C_OBJECTS = Clustal.o \
	Help.o \
	Alignment.o \
	AlignmentOutput.o \
	ObjectiveScore.o \
	Sequence.o \
	ClustalFileParser.o \
	EMBLFileParser.o \
	FileParser.o \
	FileReader.o \
	GDEFileParser.o \
	InFileStream.o \
	MSFFileParser.o \
	PearsonFileParser.o \
	PIRFileParser.o \
	RSFFileParser.o \
	ClustalWResources.o \
	DebugLog.o \
	InvalidCombination.o \
	OutputFile.o \
	Stats.o \
	SymMatrix.o \
	UserParameters.o \
	Utility.o \
	VectorOutOfRange.o \
	CommandLineParser.o \
	InteractiveMenu.o \
	Iteration.o \
	LowScoreSegProfile.o \
	MSA.o \
	MyersMillerProfileAlign.o \
	ProfileBase.o \
	ProfileStandard.o \
	ProfileWithSub.o \
	FastPairwiseAlign.o \
	FullPairwiseAlign.o \
	SubMatrix.o \
	AlignmentSteps.o \
	ClusterTree.o \
	ClusterTreeOutput.o \
	NJTree.o \
	RandomGenerator.o \
	Tree.o \
	TreeInterface.o \
	UnRootedClusterTree.o \
	Node.o \
	RootedClusterTree.o \
	RootedGuideTree.o \
	RootedTreeOutput.o \
	UPGMAAlgorithm.o

S_OBJECTS = SparseMap.o \
	Point.o \
	PointUtils.o \
	SmartKMeans.o \
	Utils.o \


mBed: Makefile mainG.cpp $(C_OBJECTS) $(S_OBJECTS) alignment/Alignment.h alignment/Sequence.h general/DebugLog.h general/Stats.h general/UserParameters.h general/Utility.h general/clustalw.h general/ClustalWResources.h fileInput/FileReader.h interface/InteractiveMenu.h interface/CommandLineParser.h sparseMap/SparseMap.h sparseMap/SmartKMeans.h 
	$(CXX) $(CXXFLAGS) -o mBed mainG.cpp $(INCLUDES) $(C_OBJECTS) $(S_OBJECTS) 



#
SparseMap.o: Makefile sparseMap/SparseMap.cpp sparseMap/SparseMap.h fileInput/FileReader.h general/SymMatrix.h tree/TreeInterface.h tree/UPGMA/RootedGuideTree.h tree/UPGMA/RootedClusterTree.h sparseMap/Point.h
	$(CXX) $(CXXFLAGS) -c sparseMap/SparseMap.cpp $(INCLUDES)

Point.o: Makefile sparseMap/Point.cpp sparseMap/Point.h
	$(CXX) $(CXXFLAGS) -c sparseMap/Point.cpp $(INCLUDES)

PointUtils.o: Makefile sparseMap/PointUtils.cpp sparseMap/PointUtils.h sparseMap/Utils.h
	$(CXX) $(CXXFLAGS) -c sparseMap/PointUtils.cpp $(INCLUDES)

SmartKMeans.o: Makefile sparseMap/SmartKMeans.cpp sparseMap/SmartKMeans.h
	$(CXX) $(CXXFLAGS) -c sparseMap/SmartKMeans.cpp $(INCLUDES)

Utils.o: Makefile sparseMap/Utils.cpp sparseMap/Utils.h sparseMap/PointUtils.h
	$(CXX) $(CXXFLAGS) -c sparseMap/Utils.cpp $(INCLUDES)


#
Clustal.o: Clustal.cpp Clustal.h pairwise/FullPairwiseAlign.h pairwise/FastPairwiseAlign.h multipleAlign/MSA.h multipleAlign/LowScoreSegProfile.h multipleAlign/Iteration.h tree/TreeInterface.h general/debuglogObject.h general/statsObject.h alignment/ObjectiveScore.h general/ClustalWResources.h Help.h
	$(CXX) $(CXXFLAGS) -c Clustal.cpp $(INCLUDES)

Help.o: Help.cpp Help.h
	$(CXX) $(CXXFLAGS) -c Help.cpp $(INCLUDES)

# alignment
Alignment.o: alignment/Alignment.cpp alignment/Alignment.h
	$(CXX) $(CXXFLAGS) -c alignment/Alignment.cpp $(INCLUDES)

AlignmentOutput.o: alignment/AlignmentOutput.cpp alignment/AlignmentOutput.h
	$(CXX) $(CXXFLAGS) -c alignment/AlignmentOutput.cpp $(INCLUDES)

ObjectiveScore.o: alignment/ObjectiveScore.cpp alignment/ObjectiveScore.h alignment/Alignment.h 
	$(CXX) $(CXXFLAGS) -c alignment/ObjectiveScore.cpp $(INCLUDES)

Sequence.o: alignment/Sequence.cpp alignment/Sequence.h
	$(CXX) $(CXXFLAGS) -c alignment/Sequence.cpp $(INCLUDES)


# fileInput
ClustalFileParser.o: fileInput/ClustalFileParser.cpp fileInput/ClustalFileParser.h 
	$(CXX) $(CXXFLAGS) -c ./fileInput/ClustalFileParser.cpp $(INCLUDES)

EMBLFileParser.o: fileInput/EMBLFileParser.cpp fileInput/EMBLFileParser.h 
	$(CXX) $(CXXFLAGS) -c fileInput/EMBLFileParser.cpp $(INCLUDES)

FileParser.o: fileInput/FileParser.cpp fileInput/FileParser.h 
	$(CXX) $(CXXFLAGS) -c fileInput/FileParser.cpp $(INCLUDES)

FileReader.o: fileInput/FileReader.cpp fileInput/FileReader.h 
	$(CXX) $(CXXFLAGS) -c fileInput/FileReader.cpp $(INCLUDES)

GDEFileParser.o: fileInput/GDEFileParser.cpp fileInput/GDEFileParser.h 
	$(CXX) $(CXXFLAGS) -c fileInput/GDEFileParser.cpp $(INCLUDES)

InFileStream.o: fileInput/InFileStream.cpp fileInput/InFileStream.h 
	$(CXX) $(CXXFLAGS) -c fileInput/InFileStream.cpp $(INCLUDES)

MSFFileParser.o: fileInput/MSFFileParser.cpp fileInput/MSFFileParser.h 
	$(CXX) $(CXXFLAGS) -c fileInput/MSFFileParser.cpp $(INCLUDES)

PearsonFileParser.o: fileInput/PearsonFileParser.cpp fileInput/PearsonFileParser.h 
	$(CXX) $(CXXFLAGS) -c fileInput/PearsonFileParser.cpp $(INCLUDES)

PIRFileParser.o: fileInput/PIRFileParser.cpp fileInput/PIRFileParser.h 
	$(CXX) $(CXXFLAGS) -c fileInput/PIRFileParser.cpp $(INCLUDES)

RSFFileParser.o: fileInput/RSFFileParser.cpp fileInput/RSFFileParser.h 
	$(CXX) $(CXXFLAGS) -c fileInput/RSFFileParser.cpp $(INCLUDES)


# general
ClustalWResources.o: general/ClustalWResources.cpp general/ClustalWResources.h
	$(CXX) $(CXXFLAGS) -c general/ClustalWResources.cpp $(INCLUDES)

DebugLog.o: general/DebugLog.cpp general/DebugLog.h
	$(CXX) $(CXXFLAGS) -c general/DebugLog.cpp $(INCLUDES)

InvalidCombination.o: general/InvalidCombination.cpp 
	$(CXX) $(CXXFLAGS) -c general/InvalidCombination.cpp $(INCLUDES)

OutputFile.o: general/OutputFile.cpp general/OutputFile.h
	$(CXX) $(CXXFLAGS) -c general/OutputFile.cpp $(INCLUDES)

Stats.o: general/Stats.cpp general/Stats.h alignment/Alignment.h 
	$(CXX) $(CXXFLAGS) -c general/Stats.cpp $(INCLUDES)

SymMatrix.o: general/SymMatrix.cpp general/SymMatrix.h
	$(CXX) $(CXXFLAGS) -c general/SymMatrix.cpp $(INCLUDES)

UserParameters.o: general/UserParameters.cpp general/UserParameters.h
	$(CXX) $(CXXFLAGS) -c general/UserParameters.cpp $(INCLUDES)

Utility.o: general/Utility.cpp general/Utility.h
	$(CXX) $(CXXFLAGS) -c general/Utility.cpp $(INCLUDES)

VectorOutOfRange.o: general/VectorOutOfRange.cpp general/VectorOutOfRange.h
	$(CXX) $(CXXFLAGS) -c general/VectorOutOfRange.cpp $(INCLUDES)

# interface
CommandLineParser.o: interface/CommandLineParser.cpp interface/CommandLineParser.h general/Utility.h
	$(CXX) $(CXXFLAGS) -c interface/CommandLineParser.cpp $(INCLUDES)

InteractiveMenu.o: interface/InteractiveMenu.cpp interface/InteractiveMenu.h
	$(CXX) $(CXXFLAGS) -c interface/InteractiveMenu.cpp $(INCLUDES)


# multipleAlign
Iteration.o: multipleAlign/Iteration.cpp multipleAlign/Iteration.h alignment/ObjectiveScore.h general/utils.h general/userparams.h tree/TreeInterface.h multipleAlign/MSA.h
	$(CXX) $(CXXFLAGS) -c multipleAlign/Iteration.cpp $(INCLUDES)

LowScoreSegProfile.o: multipleAlign/LowScoreSegProfile.cpp multipleAlign/LowScoreSegProfile.h
	$(CXX) $(CXXFLAGS) -c multipleAlign/LowScoreSegProfile.cpp $(INCLUDES)

MSA.o: multipleAlign/MSA.cpp multipleAlign/MSA.h
	$(CXX) $(CXXFLAGS) -c multipleAlign/MSA.cpp $(INCLUDES)

MyersMillerProfileAlign.o: multipleAlign/MyersMillerProfileAlign.cpp multipleAlign/MyersMillerProfileAlign.h
	$(CXX) $(CXXFLAGS) -c multipleAlign/MyersMillerProfileAlign.cpp $(INCLUDES)

ProfileBase.o: multipleAlign/ProfileBase.cpp multipleAlign/ProfileBase.h
	$(CXX) $(CXXFLAGS) -c multipleAlign/ProfileBase.cpp $(INCLUDES)

ProfileStandard.o: multipleAlign/ProfileStandard.cpp multipleAlign/ProfileStandard.h
	$(CXX) $(CXXFLAGS) -c multipleAlign/ProfileStandard.cpp $(INCLUDES)

ProfileWithSub.o: multipleAlign/ProfileWithSub.cpp multipleAlign/ProfileWithSub.h
	$(CXX) $(CXXFLAGS) -c multipleAlign/ProfileWithSub.cpp $(INCLUDES)

# pairwise
FastPairwiseAlign.o: pairwise/FastPairwiseAlign.cpp pairwise/FastPairwiseAlign.h
	$(CXX) $(CXXFLAGS) -c pairwise/FastPairwiseAlign.cpp $(INCLUDES)

FullPairwiseAlign.o: pairwise/FullPairwiseAlign.cpp pairwise/FullPairwiseAlign.h
	$(CXX) $(CXXFLAGS) -c pairwise/FullPairwiseAlign.cpp $(INCLUDES)


#substitutionMatrix
SubMatrix.o: substitutionMatrix/SubMatrix.cpp substitutionMatrix/SubMatrix.h
	$(CXX) $(CXXFLAGS) -c substitutionMatrix/SubMatrix.cpp $(INCLUDES)

# tree
AlignmentSteps.o: tree/AlignmentSteps.cpp tree/AlignmentSteps.h
	$(CXX) $(CXXFLAGS) -c tree/AlignmentSteps.cpp $(INCLUDES)

ClusterTree.o: tree/ClusterTree.cpp tree/ClusterTree.h general/OutputFile.h general/utils.h 
	$(CXX) $(CXXFLAGS) -c tree/ClusterTree.cpp $(INCLUDES)

ClusterTreeOutput.o: tree/ClusterTreeOutput.cpp tree/ClusterTreeOutput.h
	$(CXX) $(CXXFLAGS) -c tree/ClusterTreeOutput.cpp $(INCLUDES)

NJTree.o: tree/NJTree.cpp tree/NJTree.h
	$(CXX) $(CXXFLAGS) -c tree/NJTree.cpp $(INCLUDES)

RandomGenerator.o: tree/RandomGenerator.cpp tree/RandomGenerator.h
	$(CXX) $(CXXFLAGS) -c tree/RandomGenerator.cpp $(INCLUDES)

Tree.o: tree/Tree.cpp tree/Tree.h
	$(CXX) $(CXXFLAGS) -c tree/Tree.cpp $(INCLUDES)

TreeInterface.o: tree/TreeInterface.cpp tree/TreeInterface.h tree/TreeInterface-G.cpp tree/TreeInterface-G.h
	$(CXX) $(CXXFLAGS) -c tree/TreeInterface.cpp $(INCLUDES)

UnRootedClusterTree.o: tree/UnRootedClusterTree.cpp tree/UnRootedClusterTree.h general/OutputFile.h
	$(CXX) $(CXXFLAGS) -c tree/UnRootedClusterTree.cpp $(INCLUDES)

# tree/UPGMA/
Node.o: tree/UPGMA/Node.cpp tree/UPGMA/Node.h general/Utility.h
	$(CXX) $(CXXFLAGS) -c tree/UPGMA/Node.cpp $(INCLUDES)

RootedClusterTree.o: tree/UPGMA/RootedClusterTree.cpp tree/UPGMA/RootedClusterTree.h
	$(CXX) $(CXXFLAGS) -c tree/UPGMA/RootedClusterTree.cpp $(INCLUDES)

RootedGuideTree.o: tree/UPGMA/RootedGuideTree.cpp tree/UPGMA/RootedGuideTree.h
	$(CXX) $(CXXFLAGS) -c tree/UPGMA/RootedGuideTree.cpp $(INCLUDES)

RootedTreeOutput.o: tree/UPGMA/RootedTreeOutput.cpp tree/UPGMA/RootedTreeOutput.h
	$(CXX) $(CXXFLAGS) -c tree/UPGMA/RootedTreeOutput.cpp $(INCLUDES)

UPGMAAlgorithm.o: tree/UPGMA/UPGMAAlgorithm.cpp tree/UPGMA/UPGMAAlgorithm.h general/SymMatrix.h general/debuglogObject.h general/clustalw.h
	$(CXX) $(CXXFLAGS) -c tree/UPGMA/UPGMAAlgorithm.cpp $(INCLUDES)





clean:
	rm -f core clustalw2 *.o *~ lib*.a mBed

dataclean:
	rm -f *.aln *.dnd gmon.out coordinates.out  distMat.out hash.dat distPairs.dat collisions.dat


cleantree:
	make -f Makefile clean
	make -f Makefile dataclean
	cd ./tree/UPGMA;		rm -f *~; cd ../..
	cd ./tree/kmeans;		rm -f *~; cd ../..
	cd ./tree;			rm -f *~; cd ..
	cd ./pairwise;			rm -f *~; cd ..
	cd ./substitutionMatrix;	rm -f *~; cd ..
	cd ./sparseMap;			rm -f *~; cd ..
	cd ./multipleAlign;		rm -f *~; cd ..
	cd ./general;			rm -f *~; cd ..
	cd ./interface;			rm -f *~; cd ..
	cd ./fileInput;			rm -f *~; cd ..
	cd ./alignment;			rm -f *~; cd ..

