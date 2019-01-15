/*
Copyright (C) 2014 Itoh Laboratory, Tokyo Institute of Technology

This file is part of Platanus.

Platanus is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

Platanus is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with Platanus; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "scaffold.h"
#include "seqlib.h"
#include "kmer.h"
#include <cmath>
#include <algorithm>
#include <sstream>
#include <climits>
#include <cfloat>
#include <iomanip>

using std::vector;
using std::string;
using std::unordered_map;
using std::cerr;
using std::endl;


//////////////////////////////////////////////////////////////////////////////////////
// const parameter define
//////////////////////////////////////////////////////////////////////////////////////
const unsigned Scaffold::MIN_SCAFFOLD_LEN = 100;
const unsigned Scaffold::MIN_TOL_FACTOR = 2;
const unsigned Scaffold::MAX_TOL_FACTOR = 3;
const double Scaffold::DEFAULT_INS_CUTOFF_RATE = 0.05;


//////////////////////////////////////////////////////////////////////////////////////
// default constructor
//////////////////////////////////////////////////////////////////////////////////////
Scaffold::Scaffold()
: optionMinIns(), optionAveIns(), optionSDIns()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-e"] = "";
    optionSingleArgs["-s"] = "32";
    optionSingleArgs["-v"] = "32";
    optionSingleArgs["-l"] = "3";
    optionSingleArgs["-u"] = "0.1";
    optionSingleArgs["-t"] = "1";
    optionMultiArgs["-c"] = vector<string>();
    optionMultiArgs["-b"] = vector<string>();
    pairedEndSingleFileType.push_back("-ip");
    pairedEndSingleFileType.push_back("-op");
    pairedEndPairFileType.push_back("-IP");
    pairedEndPairFileType.push_back("-OP");
    optionSingleArgs["-tmp"] = ".";
}


//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void Scaffold::usage(void) const
{
    std::cerr << "\nUsage: platanus scaffold [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file (default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...]               : contig_file (fasta format)\n"
              << "    -b FILE1 [FILE2 ...]               : bubble_seq_file (fasta format)\n"
              << "    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -n{INT} INT                        : lib_id minimum_insert_size\n"
              << "    -a{INT} INT                        : lib_id average_insert_size\n"
              << "    -d{INT} INT                        : lib_id SD_insert_size\n"
              << "    -e FLOAT                           : coverage depth of homozygous region (default auto)\n"
              << "    -s INT                             : mapping seed length (default " << optionSingleArgs.at("-s") << ")\n"
              << "    -v INT                             : minimum overlap length (default " << optionSingleArgs.at("-v") << ")\n"
              << "    -l INT                             : minimum number of link (default " << optionSingleArgs.at("-l") << ")\n"
              << "    -u FLOAT                           : maximum difference for bubble crush (identity, default " << optionSingleArgs.at("-u") << ")\n"
              << "    -t INT                             : number of threads (<= " << optionSingleArgs.at("-t") << ", default 1)\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n\n\n"

              << "Outputs:\n"
              << "    PREFIX_scaffold.fa\n"
              << "    PREFIX_scaffoldBubble.fa\n"
              << "    PREFIX_scaffoldComponent.tsv\n"
             << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// initialize parameters
//////////////////////////////////////////////////////////////////////////////////////
void Scaffold::initializeParameters(void)
{
    seedLength = atoi(optionSingleArgs["-s"].c_str());
    keyLength = std::min(seedLength, platanus::ConstParam::SCAFFOLD_HASH_OVERLAP);
    bubbleThreshold = atof(optionSingleArgs["-u"].c_str());
    minLink = atoi(optionSingleArgs["-l"].c_str());
    numThread = atoi(optionSingleArgs["-t"].c_str());
    scaffoldGraph.setSeedLength(seedLength);
    scaffoldGraph.setMinOverlap(atoi(optionSingleArgs["-v"].c_str()));
    scaffoldGraph.setMinTolerenceFactor(MIN_TOL_FACTOR);

    sort(optionPairFile.begin(), optionPairFile.end());
    numFilePerLibraryID.resize(optionPairFile.size());
    libraryIDList.resize(optionPairFile.size());
    unsigned numLibrary = 0;
    for (unsigned i = 0; i < optionPairFile.size(); ++i) {
        ++(numFilePerLibraryID[numLibrary]);
        libraryIDList[numLibrary] = optionPairFile[i].libraryID;
        if (optionPairFile[i].libraryID != optionPairFile[i + 1].libraryID) {
            ++numLibrary;
        }
    }

    libraryMT.resize(numLibrary);
    omp_set_num_threads(numThread);

	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());
}


//////////////////////////////////////////////////////////////////////////////////////
// exec scaffold
//////////////////////////////////////////////////////////////////////////////////////
void Scaffold::exec(void)
{
    initializeParameters();
    mapLibraryAndInitscaffoldGraph(numThread);

	for (unsigned i = 0; i < libraryMT.size(); ++i) {
		if (i > 0 && libraryMT[i][0].getAverageInsSize() == 0) {
			vector<long> insSizeDistribution;
			libraryMT[i][0].readInsertSizeFile(insSizeDistribution);
			scaffoldGraph.insertSizeDistribution(libraryMT[i], insSizeDistribution, numThread);
			libraryMT[i][0].estimateInsSize(insSizeDistribution, libraryMT[i - 1][0].getAverageInsSize());

			std::ostringstream outStream;
			outStream << optionSingleArgs["-o"] << "_lib" << (i + 1) << "_insFreq.tsv";
			printInsertSizeFreq(outStream.str(), insSizeDistribution);
		}
		cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize()
			 << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;
		scaffoldGraph.setSeqLib(libraryMT[i]);

		for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
			scaffoldGraph.setTolerence(INT64_MAX);
			scaffoldGraph.setMinLink(std::max(minLink, scaffoldGraph.estimateLink()));
			scaffoldGraph.makeGraph(numThread);
			scaffoldGraph.setTolerence(tolerenceFactor * libraryMT[i][0].getSDInsSize());
			scaffoldGraph.removeHeteroOverlap();

			scaffoldGraph.setTolerence(tolerenceFactor * libraryMT[i][0].getSDInsSize());
			cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << scaffoldGraph.getTolerence() << endl;
			scaffoldGraph.setMinLink(std::max(minLink, scaffoldGraph.estimateLink()));
			scaffoldGraph.makeGraph(numThread);
			scaffoldGraph.crushHeteroBubble(averageCoverage);
			scaffoldGraph.crushBubbleIterative(bubbleThreshold, averageCoverage, numThread);
			scaffoldGraph.deleteErroneousEdgeIterative(numThread);
			if (i > 0) {
				scaffoldGraph.deleteRepeatEdge();
				scaffoldGraph.split();
				scaffoldGraph.makeGraph(numThread);
				scaffoldGraph.crushHeteroBubble(averageCoverage);
				scaffoldGraph.crushBubbleIterative(bubbleThreshold, averageCoverage, numThread);
				scaffoldGraph.deleteErroneousEdgeIterative(numThread);
			}
			scaffoldGraph.detectRepeat(averageCoverage);
			scaffoldGraph.makeScaffold();

			scaffoldGraph.makeGraph(numThread);
			scaffoldGraph.crushHeteroBubble(averageCoverage);
			scaffoldGraph.crushBubbleIterative(bubbleThreshold, averageCoverage, numThread);
			scaffoldGraph.deleteHeteroEdge();
			scaffoldGraph.detectRepeat(averageCoverage);
			scaffoldGraph.makeScaffold();
		}
	}

	for (unsigned i = 0; i < libraryMT.size(); ++i) {
		scaffoldGraph.setSeqLib(libraryMT[i]);
		scaffoldGraph.setMinLink(minLink);
		cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize()
			 << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;
		for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
			scaffoldGraph.setTolerence(INT64_MAX);
			scaffoldGraph.makeGraph(numThread);
			scaffoldGraph.setTolerence(tolerenceFactor * libraryMT[i][0].getSDInsSize());
			scaffoldGraph.removeHeteroOverlap();

			scaffoldGraph.setTolerence(tolerenceFactor * libraryMT[i][0].getSDInsSize());
			cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << scaffoldGraph.getTolerence() << endl;
			scaffoldGraph.makeGraph(numThread);
			scaffoldGraph.crushHeteroBubble(averageCoverage);
			scaffoldGraph.crushBubbleIterative(bubbleThreshold, averageCoverage, numThread);
			scaffoldGraph.deleteErroneousEdgeIterative(numThread);
			scaffoldGraph.deleteRepeatEdge();
			scaffoldGraph.detectRepeat(averageCoverage);
			scaffoldGraph.makeScaffold();

			scaffoldGraph.makeGraph(numThread);
			scaffoldGraph.crushHeteroBubble(averageCoverage);
			scaffoldGraph.crushBubbleIterative(bubbleThreshold, averageCoverage, numThread);
			scaffoldGraph.deleteHeteroEdge();
			scaffoldGraph.detectRepeat(averageCoverage);
			scaffoldGraph.makeScaffold();
		}
		scaffoldGraph.splitLowCoverageLinkAndDeleteErrorneousMappedPair(libraryMT, minLink, numThread);
	}

    outputAndAfterTreatment();
    cerr << "scaffold completed!" << endl;
}


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
void Scaffold::mapLibraryAndInitscaffoldGraph(const int numThread)
{
    std::unique_ptr<HeteroMapper> mapper(new HeteroMapper(seedLength, keyLength));
    platanus::Contig contig;
    platanus::Contig bubble;

    readLibrary(mapper, contig, bubble, libraryMT, numFilePerLibraryID, libraryMT.size(), numThread);
    cerr << "CONTIG_AVERAGE_COVERAGE = " << averageCoverage << endl;
    mapper->mergeBubble(numThread);

    unsigned nowFileNumber = 0;
    for (unsigned i = 0; i < libraryMT.size(); ++i) {
		libraryMT[i][0].setAverageInsSize(0);
        int nowLibraryID = optionPairFile[nowFileNumber].libraryID;
        cerr << "[LIBRARY " << libraryIDList[i] << "]" << endl;
        // set average length and minimum insert size
        // estimate insert size
        libraryMT[i][0].setAverageLength((long)((double)(libraryMT[i][0].getTotalLength()) / (2 * libraryMT[i][0].getNumPair()) + 0.5));
        long minInsertion = optionMinIns.find(nowLibraryID) == optionMinIns.end() ? 0 : optionMinIns[nowLibraryID];

        mapper->mapPairMT(libraryMT[i], minInsertion, numThread);

        libraryMT[i][0].setInsCutoffRate(optionInsCutoffRate.find(nowLibraryID) == optionInsCutoffRate.end() ? DEFAULT_INS_CUTOFF_RATE : optionInsCutoffRate[nowLibraryID]);
        if (optionAveIns.find(nowLibraryID) != optionAveIns.end() || optionSDIns.find(nowLibraryID) != optionSDIns.end()) {
            if (optionAveIns.find(nowLibraryID) != optionAveIns.end()) {
                libraryMT[i][0].setAverageInsSize(optionAveIns[nowLibraryID]);
                std::cerr << "Average insert size specified: AVE = " << libraryMT[i][0].getAverageInsSize() << std::endl;
            }
            if (optionSDIns.find(nowLibraryID) != optionSDIns.end()) {
                libraryMT[i][0].setSDInsSize(optionSDIns[nowLibraryID]);
            } else {
                libraryMT[i][0].setSDInsSize(static_cast<long>(static_cast<double>(libraryMT[i][0].getAverageInsSize()) / 10.0 + 0.5));
            }
        }
        nowFileNumber += numFilePerLibraryID[i];
    }
    scaffoldGraph.setSeqLib(libraryMT[0]);
    scaffoldGraph.initScaffolding(contig.coverage, mapper->contigMap, averageCoverage, bubbleThreshold);

    vector<long> insSizeDistribution;
    libraryMT[0][0].readInsertSizeFile(insSizeDistribution);
    scaffoldGraph.insertSizeDistribution(libraryMT[0], insSizeDistribution, numThread);

    vector<long> seqLengths;
    scaffoldGraph.scaffoldLengthList(seqLengths);

	if (libraryMT[0][0].getAverageInsSize() == 0)
		libraryMT[0][0].estimateInsSize(insSizeDistribution, 0);

    std::ostringstream outStream;
    outStream << optionSingleArgs["-o"] << "_lib" << 1 << "_insFreq.tsv";
    libraryMT[0][0].printInsertSizeFreq(outStream.str());
    cerr << "[LIBRARY " << 1 << "]\nAVE_INS = " << libraryMT[0][0].getAverageInsSize()
         << ", SD_INS = " << libraryMT[0][0].getSDInsSize() << endl;

    scaffoldGraph.saveOverlap(mapper->contigMap, platanus::ConstParam::SCAFFOLD_HASH_OVERLAP, libraryMT[0][0].getSDInsSize() * MIN_TOL_FACTOR, numThread);
    scaffoldGraph.countBubble(bubble, mapper->getBubblePositionRef());
    scaffoldGraph.classifyNode();
    cerr << "destructing mapper objects..." << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////
// Read library files
// Change read function which fasta or fastq, single file or pair one.
//////////////////////////////////////////////////////////////////////////////////////
void Scaffold::readLibrary(std::unique_ptr<HeteroMapper> &mapper, platanus::Contig &contig, platanus::Contig &bubble, vector<vector<SeqLib> > &libraryMT, vector<int> &numFilePerLibraryID, const int numLibrary, const int numThread)
{
    omp_set_num_threads(numThread);
    # pragma omp parallel for schedule(static, 1)
    for (int i = -2; i < numLibrary; ++i) {
        try {
            // load contig file
            if (i == -2) {
                for (unsigned i = 0; i < optionMultiArgs["-c"].size(); ++i) {
                    contig.readFastaCoverage(optionMultiArgs["-c"][i]);
                }

				if (optionSingleArgs["-e"] == "")
					averageCoverage = contig.calculateAverageCoverageExcludingOutlier(contig.getSeqLengthMedian());
				else
					averageCoverage = atof(optionSingleArgs["-e"].c_str());

                mapper->setContigMap(contig);
                mapper->makeKmerTableContigMap();
            // load contigBubble file
            } else if (i == -1) {
                for (unsigned i = 0; i < optionMultiArgs["-b"].size(); ++i) {
                    bubble.readFastaCoverage(optionMultiArgs["-b"][i]);
                }
                mapper->setBubbleMap(bubble);
                mapper->makeKmerTableBubbleMap();
            // load read file
            } else {
                unsigned nowFileNumber = 0;
                libraryMT[i].resize(numThread);
                for (int j = 0; j < i; ++j)
                    nowFileNumber += numFilePerLibraryID[j];

                for (int j = 0; j < numThread; ++j) {
                    libraryMT[i][j].makeTempPairFP();
                }
                for (int j = 0; j < numFilePerLibraryID[i]; ++j) {
                    bool isFastq, isMate;
                    isMate = optionPairFile[nowFileNumber+j].libraryType == "-op"
                           || optionPairFile[nowFileNumber+j].libraryType == "-OP" ? true : false;
                    // if read file is single file (forward and reverse seq contain same file)
                    if (optionPairFile[nowFileNumber+j].libraryType == "-ip" || optionPairFile[nowFileNumber+j].libraryType == "-op") {
                        platanus::FILETYPE fileFormat = checkFileFormat(optionPairFile[nowFileNumber+j].fileFirst);
                        switch (fileFormat) {
                        case platanus::UNKNOWN:
                            throw platanus::FormatError();
                            break;
                        case platanus::FASTA:
                            isFastq = false;
                            break;
                        default:
                            isFastq = true;
                            break;
                        }
                        ReadFastaSingleMT(libraryMT[i], optionPairFile[nowFileNumber+j].fileFirst, numThread, isMate, isFastq);
                    // if read file is pair file (forward and reverse seq contain differenct file)
                    } else {
                        platanus::FILETYPE fileFormat1 = checkFileFormat(optionPairFile[nowFileNumber+j].fileFirst);
                        platanus::FILETYPE fileFormat2 = checkFileFormat(optionPairFile[nowFileNumber+j].fileSecond);
                        if (fileFormat1 != fileFormat2) {
                            throw platanus::FormatError("Different file type in paired-file.");
                        }
                        switch (fileFormat1) {
                        case platanus::UNKNOWN:
                            throw platanus::FormatError();
                            break;
                        case platanus::FASTA:
                            isFastq = false;
                            break;
                        default:
                            isFastq = true;
                            break;
                        }
                        ReadFastaPairMT(libraryMT[i], optionPairFile[nowFileNumber+j].fileFirst, optionPairFile[nowFileNumber+j].fileSecond, numThread, isMate, isFastq);
                    }
                }
            }
        } catch (platanus::ErrorBase &e) {
            e.showErrorMessage();
            exit(e.getID());
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// destroy graph
//////////////////////////////////////////////////////////////////////////////////////
void Scaffold::outputAndAfterTreatment(void)
{
    cerr << "writing scaffold files..." << endl;
    string outFilename(optionSingleArgs["-o"]);
    string componentFilename(outFilename);
    outFilename += "_scaffold.fa";
    componentFilename += "_scaffoldComponent.tsv";
    scaffoldGraph.cutAndPrintSeq(MIN_SCAFFOLD_LEN, outFilename, componentFilename);

    outFilename = optionSingleArgs["-o"] + "_scaffoldBubble.fa";
    scaffoldGraph.printScaffoldBubble(outFilename);
}



void Scaffold::printInsertSizeFreq(const std::string &outputFilename, const vector<long>& distribution)
{
    long maxInsSize = distribution.size() - 1;
    long i = 1;
    std::ofstream out(outputFilename.c_str());

    for (; i <= maxInsSize; ++i) {
        out << i << "\t" << distribution[i] << std::endl;
    }
    unsigned long long tmp = 0;
    for (; static_cast<unsigned long long>(i) < distribution.size(); ++i) {
        tmp += distribution[i];
    }
    out << maxInsSize + 1 << "\t" << tmp << std::endl;

    out.close();
}
