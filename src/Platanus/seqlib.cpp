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

#include "seqlib.h"
#include <vector>
#include <cmath>

using std::vector;
using std::cerr;
using std::endl;
using std::ifstream;


//////////////////////////////////////////////////////////////////////////////////////
// static constant
//////////////////////////////////////////////////////////////////////////////////////
//const double SeqLib::INS_DISTR_TRUNC = 0.01;
const double SeqLib::INS_DISTR_TRUNC = 0.025;
const double SeqLib::INS_DISTR_TRUNC_SD_RATE = 3.0;
const unsigned SeqLib::INS_DISTR_TRUNC_NUM_ITERATION = 1000;
const double SeqLib::INS_CUTOFF_RATE_TO_PEAK = 0.5;
const long SeqLib::INS_PEAK_WINDOW = 101;



//////////////////////////////////////////////////////////////////////////////////////
// cut outlier value in distribution of insert size
// default cut fraction is 0.01 in both-side
//////////////////////////////////////////////////////////////////////////////////////
template <typename T>
void SeqLib::truncateDistribution(vector<T> &distribution, const double edge) const
{
    if (distribution.size() <= 1)
        return;

    T totalSum = 0;
    for (unsigned long i = 1; i < distribution.size(); ++i)
        totalSum += distribution[i] * i;

    T edgeSum = 0;
    for (unsigned long i = 1; static_cast<double>(edgeSum) < totalSum * edge; ++i) {
        edgeSum += distribution[i] * i;
        distribution[i - 1] = 0;
    }
    edgeSum = (distribution.size() - 1) * distribution[distribution.size() - 1];
    for (unsigned long i = distribution.size() - 2; static_cast<double>(edgeSum) < totalSum * edge; --i) {
        edgeSum += distribution[i] * i;
        distribution[i + 1] = 0;
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// cut outlier value in distribution of insert size
// default cut fraction is 0.01 in both-side
//////////////////////////////////////////////////////////////////////////////////////
template <typename T>
void SeqLib::truncateDistributionByNumber(vector<T> &distribution, const double edge) const
{
    if (distribution.size() <= 1)
        return;

    T totalSum = 0;
    for (unsigned long i = 1; i < distribution.size(); ++i)
        totalSum += distribution[i];
    T finalSum = edge / 2.0 * totalSum;

    unsigned long lowerThreshold = 0;
    T edgeSum = 0;
    for (; edgeSum < finalSum; ++lowerThreshold) {
        edgeSum += distribution[lowerThreshold];
        if (edgeSum > finalSum) {
            distribution[lowerThreshold] = edgeSum - finalSum;
            break;
        } else {
            distribution[lowerThreshold] = 0;
        }
    }
    unsigned long upperThreshold = distribution.size() - 1;
    edgeSum = 0;
    for (; edgeSum < finalSum; --upperThreshold) {
        edgeSum += distribution[upperThreshold];
        if (edgeSum + distribution[upperThreshold] > finalSum) {
            distribution[upperThreshold] = edgeSum - finalSum;
            break;
        } else {
            distribution[upperThreshold] = 0;
        }
    }
    cerr << "LOWER_THRESHOLD = " << lowerThreshold << ", UPPER_THRESHOLD = " << upperThreshold << endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// cut outlier value in distribution of insert size
//////////////////////////////////////////////////////////////////////////////////////
template <typename T>
void SeqLib::truncateDistributionBySD(vector<T> &distribution, const double edge)
{
    if (distribution.size() <= 1)
        return;

    averageInsSize = calcDistributionAverage(distribution, distribution.size() - 1) + 0.5;
    sdInsSize = calcDistributionSD(distribution, distribution.size() - 1) + 0.5;
    long bestAverage = getAverageInsSize();
    long bestSD = getSDInsSize();

    long lowerThreshold;
    long upperThreshold;
    unsigned count = 0;
    while (true) {
        long preAverage = getAverageInsSize();
        long preSD = getSDInsSize();
        vector<T> tmpDistribution = distribution;
        lowerThreshold = preAverage - (edge * preSD - 0.5);
        upperThreshold = preAverage + (edge * preSD + 0.5);
        for (long idx = 0; idx < lowerThreshold; ++idx) {
            tmpDistribution[idx] = 0;
        }
        for (unsigned long idx = upperThreshold; idx < tmpDistribution.size(); ++idx) {
            tmpDistribution[idx] = 0;
        }

        averageInsSize = calcDistributionAverage(distribution, distribution.size() - 1) + 0.5;
        sdInsSize = calcDistributionSD(distribution, distribution.size() - 1) + 0.5;
        long tmpAverage = getAverageInsSize();
        long tmpSD = getSDInsSize();
        if (tmpAverage == preAverage) {
            break;
        }
        if (tmpSD <= preSD) {
            bestAverage = tmpAverage;
            bestSD = tmpSD;
        }
        if (++count > INS_DISTR_TRUNC_NUM_ITERATION) {
	  cerr << "Warning: fail to estimate AVE_INS and SD_INS!" << endl;
	  lowerThreshold = bestAverage - (edge * bestSD + 0.5);
	  upperThreshold = bestAverage + (edge * bestSD + 0.5);
	  break;
        }
    }
    cerr << "LOWER_THRESHOLD = " << lowerThreshold << ", UPPER_THRESHOLD = " << upperThreshold << endl;

    for (long idx = 0; idx < lowerThreshold; ++idx) {
        distribution[idx] = 0;
    }
    for (unsigned long idx = upperThreshold; idx < distribution.size(); ++idx) {
        distribution[idx] = 0;
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// calculate average of distribution
//////////////////////////////////////////////////////////////////////////////////////
template <typename T>
double SeqLib::calcDistributionAverage(const vector<T> &distribution, const long lowerBound, const long upperBound) const
{
    T sumDistribution = 0;
    T numDistribution = 0;

    for (long i = lowerBound; i <= upperBound; ++i) {
        sumDistribution += distribution[i] * i;
        numDistribution += distribution[i];
    }
    return static_cast<double>(sumDistribution) / numDistribution + 0.5;
}


//////////////////////////////////////////////////////////////////////////////////////
// calculate standard deviation of distribution
//////////////////////////////////////////////////////////////////////////////////////
template <typename T>
double SeqLib::calcDistributionSD(const vector<T> &distribution, const long lowerBound, const long upperBound) const
{
    double numerator = 0;
    T denominator = 0;

    for (long i = lowerBound; i <= upperBound; ++i) {
        numerator += (static_cast<double>(i) - averageInsSize) * (static_cast<double>(i) - averageInsSize) * distribution[i];
        denominator += distribution[i];
    }
    if (denominator > 1)
        return sqrt(numerator / static_cast<double>(denominator - 1));
    else
        return 0.0;
}


void SeqLib::normalizeDistribution(const vector<long>& preDist, vector<double>& postDist, const vector<long>& seqLengths) const
{
    vector<long> probability(preDist.size(), 0);
    for (unsigned long idx = 0; idx < seqLengths.size(); ++idx) {
        unsigned long tmpLength = seqLengths[idx];
        for (unsigned long i = 1, end = std::min(tmpLength + 1, preDist.size()); i < end; ++i) {
            probability[i] += tmpLength + 1 - i;
        }
    }

    postDist.clear();
    postDist.resize(preDist.size(), 0);
    for (unsigned long idx = 1; idx < preDist.size(); ++idx) {
        postDist[idx] = static_cast<double>(preDist[idx]) / probability[idx];
    }

    double sumPreDist  = 0;
    double sumPostDist = 0;
    for (unsigned long idx = 1; idx < preDist.size(); ++idx) {
        sumPreDist  += preDist[idx];
        sumPostDist += postDist[idx];
    }
    double rate = static_cast<double>(sumPreDist) / sumPostDist;
    for (unsigned long idx = 1; idx < postDist.size(); ++idx) {
        postDist[idx] *= rate;
    }

}


//////////////////////////////////////////////////////////////////////////////////////
// estimate library insert size
//////////////////////////////////////////////////////////////////////////////////////
void SeqLib::estimateInsSize(const vector<long>& distribution, const long minPeakThreshold)
{
	const double upperBoundFactor = 1.75;
	const double lowerBoundFactor = 0.25;

    cerr << "estimating insert-size..." << endl;

	long peakInsSize = findDistributionPeak(distribution, INS_PEAK_WINDOW, minPeakThreshold);
	long upperBound = std::min(static_cast<long>(upperBoundFactor * peakInsSize + 0.5), static_cast<long>(distribution.size() - 1));
	long lowerBound = static_cast<long>(lowerBoundFactor * peakInsSize + 0.5);

    averageInsSize = calcDistributionAverage(distribution, lowerBound, upperBound) + 0.5;
    sdInsSize = calcDistributionSD(distribution, lowerBound, upperBound) + 0.5;

    cerr << "PEAK = " << peakInsSize << endl;
    cerr << "LOWER_LIMIT (permissible range to estimate AVE_INS)= " << lowerBound << endl;
    cerr << "UPPER_LIMIT (permissible range to estimate AVE_INS)= " << upperBound << endl;
    cerr << "AVE_INS = " << averageInsSize << endl;
    cerr << "SD_INS = " << sdInsSize << endl;
}

//////////////////////////////////////////////////////////////////////////////////////
// estimate library insert size
//////////////////////////////////////////////////////////////////////////////////////
void SeqLib::estimateInsSizeNormalized(const vector<long>& insSizeDistribution, vector<long>& seqLengths)
{
    vector<double> distribution;
    normalizeDistribution(insSizeDistribution, distribution, seqLengths);
    truncateDistributionByNumber(distribution, insCutoffRate);

    //long maxInsSize = std::min(2 * findDistributionPeak(normalizedDistribution, INS_PEAK_WINDOW), static_cast<long>(normalizedDistribution.size() - 1), 0);
    averageInsSize = calcDistributionAverage(distribution, 0, distribution.size() - 1) + 0.5;
    sdInsSize = calcDistributionSD(distribution, 0, distribution.size() - 1) + 0.5;
}


void SeqLib::readInsertSizeFile(vector<long>& distribution)
{
    distribution.clear();
    long insLen;
    rewind(insertLengthFP);
    while (fread(&insLen, sizeof(long), 1, insertLengthFP)) {
        if (insLen >= static_cast<long>(distribution.size())) {
            distribution.resize(insLen + 1, static_cast<long>(0));
        }
        ++distribution[insLen];
    }
    if (distribution.size() == 0) {
        throw platanus::MapError("No read mapped in the same contig!!");
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// Output insert size frequency file
//////////////////////////////////////////////////////////////////////////////////////
void SeqLib::printInsertSizeFreq(const std::string &outputFilename)
{
    vector<long> distribution;
    readInsertSizeFile(distribution);

    std::ofstream out(outputFilename.c_str());
//    long maxInsSize = std::min(2 * findDistributionPeak(distribution, INS_PEAK_WINDOW), static_cast<long>(distribution.size() - 1));
    long maxInsSize = distribution.size() - 1;
    long i = 1;
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


//////////////////////////////////////////////////////////////////////////////////////
// find distribution peak value
//////////////////////////////////////////////////////////////////////////////////////
template <typename T>
long SeqLib::findDistributionPeak(const vector<T> &distribution, const long windowSize, long minPeakThreshold) const
{
    if (distribution.size() <= static_cast<unsigned long long>(std::min(windowSize, minPeakThreshold))) {
        return distribution.size() / 2;
    }

	minPeakThreshold = std::max(minPeakThreshold, windowSize / 2);

    T pre = 0;
    for (long i = minPeakThreshold - windowSize / 2; i < windowSize; ++i) {
        pre += distribution[i];
    }

    T peak = pre;
    long peakI = minPeakThreshold;

    for (long i = minPeakThreshold - windowSize / 2 + 1, n = distribution.size() - windowSize; i <= n; ++i) {
        T cur = pre - distribution[i - 1] + distribution[i + windowSize - 1];
        if (cur > peak) {
            peak = cur;
            peakI = i + windowSize / 2;
        }
        pre = cur;
    }

    return peakI;
}


//////////////////////////////////////////////////////////////////////////////////////
// Read single FASTA(Q) file and write tempfile
//////////////////////////////////////////////////////////////////////////////////////
void ReadFastaSingleMT(vector<SeqLib> &lib, const std::string &filename, const int numThread, const bool isMate, const bool isFastq, const bool notPair)
{
    long i;
    std::ifstream forwardFile(filename.c_str());
    if (!forwardFile) throw platanus::FILEError(filename);

    platanus::SEQ forwardSeq, reverseSeq;
    void (* const ReadFunction)(platanus::SEQ &seq, std::ifstream &ifs) = isFastq ? ReadFastqSeq : ReadFastaSeq;
    void (* const ReadHeaderFunction)(std::ifstream &ifs) = isFastq ? ReadFastqHeader : ReadFastaHeader;
    for (i = 0; i < numThread; ++i)
        fseek(lib[i].pairFP, 0, SEEK_END);
    if (!forwardFile.good()) {
        throw platanus::IOError();
    }
    ReadHeaderFunction(forwardFile);
    i = 0;
    while (forwardFile && !forwardFile.eof()) {
        ReadFunction(forwardSeq, forwardFile);
        if (!notPair && forwardFile.eof()) {
            throw platanus::FormatError("the number of read is odd in file.");
        }
        ReadFunction(reverseSeq, forwardFile);
        if (isMate) {
            forwardSeq.reverse();
            reverseSeq.reverse();
        }
        forwardSeq.writeTemporaryFile(lib[i].pairFP);
        reverseSeq.writeTemporaryFile(lib[i].pairFP);
        i = (i + 1) % numThread;
        lib[0].addTotalLength(forwardSeq.length + reverseSeq.length);
        lib[0].addNumPair(1);
    }


    forwardFile.close();

}



//////////////////////////////////////////////////////////////////////////////////////
// Read pair FASTA(Q) file and write tempfile
//////////////////////////////////////////////////////////////////////////////////////
void ReadFastaPairMT(vector<SeqLib> &lib, const std::string &filename1, const std::string &filename2, const int numThread, const bool isMate, const bool isFastq)
{
    long i;
    ifstream forwardFile(filename1.c_str());
    if (!forwardFile) throw platanus::FILEError(filename1);
    ifstream reverseFile(filename2.c_str());
    if (!reverseFile) throw platanus::FILEError(filename2);
    platanus::SEQ forwardSeq, reverseSeq;
    void (*ReadFunction)(platanus::SEQ &seq, ifstream &ifs) = isFastq ? ReadFastqSeq : ReadFastaSeq;
    void (*ReadHeaderFunction)(ifstream &ifs) = isFastq ? ReadFastqHeader : ReadFastaHeader;
    for (i = 0; i < numThread; ++i)
        fseek(lib[i].pairFP, 0, SEEK_END);
    if (!forwardFile.good() || !reverseFile.good()) {
        throw platanus::IOError();
    }
    ReadHeaderFunction(forwardFile);
    ReadHeaderFunction(reverseFile);

    i = 0;
    while (forwardFile && reverseFile && !forwardFile.eof() && !reverseFile.eof()) {
        ReadFunction(forwardSeq, forwardFile);
        ReadFunction(reverseSeq, reverseFile);
        if (isMate) {
            forwardSeq.reverse();
            reverseSeq.reverse();
        }
        forwardSeq.writeTemporaryFile(lib[i].pairFP);
        reverseSeq.writeTemporaryFile(lib[i].pairFP);
        i = (i + 1) % numThread;

        lib[0].addTotalLength(forwardSeq.length + reverseSeq.length);
        lib[0].addNumPair(1);
    }

    if (!forwardFile.eof() || !reverseFile.eof()) {
        throw platanus::FormatError("the number of read is different in paired-file.");
    }

    forwardFile.close();
    reverseFile.close();

}



//////////////////////////////////////////////////////////////////////////////////////
// Jump fasta head header
//////////////////////////////////////////////////////////////////////////////////////
void ReadFastaHeader(ifstream &ifs)
{
    std::string line;
    while (ifs && getline(ifs, line)) {
        if (line[0] == '>') {
            break;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////
// Jump fastq head header
//////////////////////////////////////////////////////////////////////////////////////
void ReadFastqHeader(ifstream &ifs)
{
    std::string line;
    while (ifs && getline(ifs, line)) {
        if (line[0] == '@') {
            break;
        }
    }
}



//////////////////////////////////////////////////////////////////////////////////////
// Read one fasta seq
//////////////////////////////////////////////////////////////////////////////////////
void ReadFastaSeq(platanus::SEQ &seq, std::ifstream &ifs)
{
    std::string line;
    seq.length = 0;
    seq.base = "";
    while (ifs && getline(ifs, line)) {
        if (line[0] == '>') {
            break;
        } else {
            seq.put(line);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////
// Read one fastq seq
//////////////////////////////////////////////////////////////////////////////////////
void ReadFastqSeq(platanus::SEQ &seq, std::ifstream &ifs)
{
    unsigned readLine = 0;
    std::string line;
    seq.length = 0;
    seq.base = "";
    while (ifs && getline(ifs, line)) {
        if (line[0] == '+') {
            break;
        } else {
            ++readLine;
            seq.put(line);
        }
    }
    for (unsigned i = 0; i < readLine + 1; ++i)
        getline(ifs, line);
}


void SeqLib::deleteErrorneousMappedPair(const std::unordered_map<std::pair<int, int>, bool, platanus::PairHash, platanus::PairEqual> &errorLink)
{
    platanus::Position forwardResult;
    platanus::Position reverseResult;
    FILE *newMappedFP = platanus::makeTemporaryFile();

    rewind(mappedFP);
    while (fread(&forwardResult, sizeof(platanus::Position), 1, mappedFP)) {
        fread(&reverseResult, sizeof(platanus::Position), 1, mappedFP);

        std::pair<int, int> key;
        if (abs(forwardResult.id) < abs(reverseResult.id)) {
            key.first = forwardResult.id;
            key.second = -(reverseResult.id);
        }
        else {
            key.first = reverseResult.id;
            key.second = -(forwardResult.id);
        }

        if (errorLink.find(key) == errorLink.end()) {
            fwrite(&forwardResult, sizeof(platanus::Position), 1, newMappedFP);
            fwrite(&reverseResult, sizeof(platanus::Position), 1, newMappedFP);
        }
    }

    fclose(mappedFP);
    mappedFP = newMappedFP;
}



