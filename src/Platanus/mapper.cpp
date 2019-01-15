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

#include "mapper.h"
#include <omp.h>

using std::vector;
using std::cerr;
using std::endl;


//////////////////////////////////////////////////////////////////////////////////////
// initilize mapper (not constructor)
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::setContig(platanus::Contig &contig)
{
    numSeq = contig.numSeq;
    seq = std::move(contig.seq);
    for (auto it = seq.begin(), end = seq.end(); it != end; ++it)
        seqPoolSize += it->base.size();
}





//////////////////////////////////////////////////////////////////////////////////////
// insert contigs into kmer hash table
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::insertHashKey(void)
{
    unsigned long long seqLength = 0;
    unsigned long long hashKey = 0;
    unsigned long long start = 0;
    unsigned long long total = 0;
    bool keyInit = false;
    unsigned j = 0;
    FILE *tmpFP = platanus::makeTemporaryFile();
    // make this struct
    platanus::Position position;
    table.clear();
    
    for (int i = 0; i < numSeq; ++i) {
        hashKey = 0;
        seqLength = seq[i].length;
        if (static_cast<long>(seqLength) < keyLength) continue;
        start = 0;
        keyInit = true;
        while (start < seqLength - keyLength + 1) {
            // initialize hash key
            if (keyInit) {
                hashKey = 0;
                for (j = 0; j < static_cast<unsigned>(keyLength - 1); ++j) {
                    // if sequence base contain "N", don't make hash key
                    if (seq[i].base[start + j] == 4) break;
                    hashKey = (hashKey << 2) | static_cast<unsigned long long>(seq[i].base[start + j]);
                }
                if (j == static_cast<unsigned>(keyLength - 1)) {
                    keyInit = false;
                } else {
                    start += j + 1;
                    continue;
                }
            }
            if (seq[i].base[start + keyLength - 1] == 4) {
                start += keyLength;
                keyInit = true;
                continue;
            }

            hashKey = ((hashKey << 2) & mask) | static_cast<unsigned long long>(seq[i].base[start + keyLength - 1]);
            position.id = i + 1;
            position.offset = start;
		//this line add key to table,not but .position
            ++(table[hashKey].num);
            fwrite(&hashKey, sizeof(unsigned long long), 1, tmpFP);
            fwrite(&position, sizeof(platanus::Position), 1, tmpFP);
            ++total;
            ++start;
        }
    }
    platanus::Position *pos = reservePositionPool(total);

    fillPositionPool(tmpFP, pos);
}



//////////////////////////////////////////////////////////////////////////////////////
// new positionPool using smart pointer
//////////////////////////////////////////////////////////////////////////////////////
platanus::Position *Mapper::reservePositionPool(const unsigned long long total)
{
    positionPool.reset(new platanus::Position[total]);
    return positionPool.get();
}



//////////////////////////////////////////////////////////////////////////////////////
// fill positionPool and complete map table
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::fillPositionPool(FILE *tmpFP, platanus::Position *pos)
{
    unsigned long long hashKey;
    platanus::Position position;
        rewind(tmpFP);
        while (fread(&hashKey, sizeof(unsigned long long), 1, tmpFP)) {
            fread(&position, sizeof(platanus::Position), 1, tmpFP);
            auto it = table.find(hashKey);
            if (it->second.position == NULL) {
                it->second.position = pos;
                pos += it->second.num;
                it->second.num = 1;
                *(it->second.position) = position;
            } else {
                *(it->second.position + it->second.num) = position;
                ++(it->second.num);
            }
        }
        fclose(tmpFP);
}



//////////////////////////////////////////////////////////////////////////////////////
// make kmer table and calculate max occurrence value
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::makeKmerTable(void)
{
    cerr << "K=" << keyLength << ", making hash table..." << endl;
    this->insertHashKey();
}



//////////////////////////////////////////////////////////////////////////////////////
// mapping kmer (however, only k <= 32)
//////////////////////////////////////////////////////////////////////////////////////
long Mapper::mapSeed(const Kmer31 &kmer, vector<platanus::Position> &positionBuffer) const
{
    positionBuffer.clear();
    // map forward kmer in table
    auto tableIt = table.find(kmer.forward);
    if (tableIt != table.end()) {
        positionBuffer.resize(tableIt->second.num);
        std::copy(tableIt->second.position, tableIt->second.position + tableIt->second.num, positionBuffer.begin());
    }

    unsigned forwardSize = positionBuffer.size();
    // map reverse kmer in table
    tableIt = table.find(kmer.reverse);
    if (tableIt != table.end()) {
        positionBuffer.insert(positionBuffer.end(), tableIt->second.position, tableIt->second.position + tableIt->second.num);
    }

    // in mapped reverse kmer, contig id *= -1
    for (unsigned i = forwardSize; i < positionBuffer.size(); ++i)
        positionBuffer[i].id *= -1;
    return positionBuffer.size();

}

long Mapper::mapKey(const Kmer31 &kmer, vector<platanus::Position> &positionBuffer) const
{ 
    positionBuffer.clear();
    unsigned int MultiCutoff = 10;
    auto tableIt = table.find(kmer.forward);
    if (tableIt != table.end()) {
	  if (tableIt->second.num > MultiCutoff-1 && seedLength==keyLength)
	  {
		positionBuffer.resize(MultiCutoff);
		std::copy(tableIt->second.position, tableIt->second.position + MultiCutoff, positionBuffer.begin());
	  }
	  else{
	  positionBuffer.resize(tableIt->second.num);
	  std::copy(tableIt->second.position, tableIt->second.position + tableIt->second.num, positionBuffer.begin());
    }
    }
    unsigned forwardSize = positionBuffer.size();
    tableIt = table.find(kmer.reverse);
    if (tableIt != table.end()) {
	  if (tableIt->second.num > MultiCutoff-1 && seedLength==keyLength)
		positionBuffer.insert(positionBuffer.end(), tableIt->second.position, tableIt->second.position + MultiCutoff);
	  else
	  positionBuffer.insert(positionBuffer.end(), tableIt->second.position, tableIt->second.position + tableIt->second.num);
    }

    for (unsigned i = forwardSize; i < positionBuffer.size(); ++i)
	  positionBuffer[i].id *= -1;
    return positionBuffer.size();

}



//////////////////////////////////////////////////////////////////////////////////////
// Merge bubble and contig by mapping
//////////////////////////////////////////////////////////////////////////////////////
void HeteroMapper::mergeBubble(const long numThread)
{
    long threadID = 0;
    long seqID = 0;
    long bufferSize = 0;
    long maxLeftLength = 0;
    long maxRightLength = 0;
    vector<platanus::Position> positionBuffer;
    platanus::Position leftResult;
    platanus::Position rightResult;
    if (bubbleMap.getNumSeq() == 0) return;

    cerr << "mapping bubbles on contigs..." << endl;
    if (bubblePosition.size() != 0)
        bubblePosition.clear();

    bubblePosition.resize(bubbleMap.getNumSeq());
    omp_set_num_threads(numThread);

    #    pragma omp parallel for  schedule(static, 1) \
        private(seqID, positionBuffer, bufferSize, maxLeftLength, maxRightLength) \
        firstprivate(leftResult, rightResult)

    for (threadID = 0; threadID < numThread; ++threadID) {
        long bubbleNumSeq = bubbleMap.getNumSeq();
        for (seqID = threadID; seqID < bubbleNumSeq; seqID += numThread) {
            Kmer31 leftKmer(keyLength);
            Kmer31 rightKmer(keyLength);
            const platanus::SEQ &bubbleSeq = bubbleMap.getSeq(seqID);
            if (bubbleSeq.length < 2 * keyLength) continue;

            if (leftKmer.setKmer(bubbleSeq, 0, keyLength, seedLength) || rightKmer.setKmer(bubbleSeq, bubbleSeq.length - keyLength, keyLength, seedLength)) continue;
            // mapping bubble left seed in contig
            bufferSize = contigMap.mapSeed(leftKmer, positionBuffer);
            maxLeftLength = 0;
            for (long i = 0; i < bufferSize; ++i) {
                long j;
                // position.id > 0 means seed is mapped forward contig
                if (positionBuffer[i].id > 0) {
                    for (j = keyLength; j < bubbleSeq.length; ++j) {
                        if (positionBuffer[i].offset + j >= contigMap.getSeqLength(positionBuffer[i].id - 1)
                            || contigMap.getBase(positionBuffer[i].id - 1, positionBuffer[i].offset + j) != bubbleSeq.base[j])
                            break;
                    }
                }
                // position.id < 0 means seed is mapped reverse contig
                else {
                    for (j = keyLength; j < bubbleSeq.length; ++j) {
                        if (positionBuffer[i].offset + keyLength - j <= 0  
                            || contigMap.getBase(-(positionBuffer[i].id) - 1, positionBuffer[i].offset + keyLength - j - 1) != (0x3^bubbleSeq.base[j]))
                            break;
                    }
                    positionBuffer[i].offset += keyLength - 1;
                }
                if (j > maxLeftLength) {
                    maxLeftLength = j;
                    leftResult =  positionBuffer[i];
                } else if (j == maxLeftLength)
                    leftResult.id = 0;
            }

            if (maxLeftLength < seedLength) continue;

            // mapping bubble right seed in contig
            bufferSize = contigMap.mapSeed(rightKmer, positionBuffer);
            maxRightLength = 0;
            for (long i = 0; i < bufferSize; ++i) {
                long j;
                // position.id > 0 means seed is mapped forward contig
                if (positionBuffer[i].id > 0) {
                    for (j = keyLength; j < bubbleSeq.length; ++j) {
                        if (positionBuffer[i].offset + keyLength - j <= 0
                            || contigMap.getBase(positionBuffer[i].id - 1, positionBuffer[i].offset + keyLength - j - 1) != bubbleSeq.base[bubbleSeq.length - j - 1])
                            break;
                    }
                    positionBuffer[i].offset += keyLength - 1;
                }
                // position.id < 0 means seed is mapped reverse contig
                else {
                    for (j = keyLength; j < bubbleSeq.length; ++j) {
                        if (positionBuffer[i].offset + j >= contigMap.getSeqLength(-(positionBuffer[i].id) - 1)
                            || contigMap.getBase(-(positionBuffer[i].id) - 1, positionBuffer[i].offset + j) != (0x3^bubbleSeq.base[bubbleSeq.length - j - 1]))
                            break;
                    }
                }
                if (j > maxRightLength) {
                    maxRightLength = j;
                    rightResult =  positionBuffer[i];
                } else if (j == maxRightLength)
                    rightResult.id = 0;
            }

            if (maxRightLength < seedLength || rightResult.id != leftResult.id || leftResult.id == 0 || rightResult.id == 0) continue;

            bubblePosition[seqID].id = leftResult.id;
            bubblePosition[seqID].offset = (leftResult.offset + rightResult.offset) / 2;
        }
    }
}




//////////////////////////////////////////////////////////////////////////////////////
// decide the contig where the variable read is mapped
// this version, I only translate C language.
// it can be more acceleratable using kmer mapping over k > 32
//////////////////////////////////////////////////////////////////////////////////////
platanus::Position Mapper::mapRead(const platanus::SEQ &read, vector<platanus::Position> &positionBuffer)
{
    long numResult = 0;
    // vector<platanus::Position> result(platanus::ConstParam::MAX_READ_LEN);
    vector<platanus::Position> result(read.length/seedLength + 2);
    Kmer31 kmer;
    long bufferSize;
    long firstMax, secondMax;
    int i, j, k = 0;
    for (i = read.length - seedLength; i > -(seedLength); i -= seedLength) {
        if (i < 0)
            i = 0;
        kmer.forward = kmer.reverse = 0;
        for (j = 0; j < keyLength; ++j) {
            if (read.base[i + j] == 4) break;
            kmer.forward = (kmer.forward << 2) | static_cast<unsigned long long>(read.base[i + j]);
            kmer.reverse = (kmer.reverse << 2) | static_cast<unsigned long long>(0x3^read.base[i + keyLength - j -1]);
        }

        if (j != keyLength) continue;

        // first, mapping kmer <= 32
        bufferSize = mapKey(kmer, positionBuffer);
        result[numResult].id = 0;
        // extend mapping kmer >= 32
        for (j = 0; j < bufferSize; ++j) {
            // position.id > 0 means seed is mapped forward contig
            if (positionBuffer[j].id > 0) {
                if (positionBuffer[j].offset > seq[positionBuffer[j].id - 1].length - seedLength)
                    continue;
                for (k = keyLength; k < seedLength; ++k) {
                    if (seq[positionBuffer[j].id - 1].base[positionBuffer[j].offset + k] != read.base[i + k])
                        break;
                }
                positionBuffer[j].offset -= i;
            }
            // position.id < 0 means seed is mapped reverse contig
            else {
                if (positionBuffer[j].offset < seedLength - keyLength)
                    continue;
                for (k = keyLength; k < seedLength; ++k) {
                    if (seq[-(positionBuffer[j].id) - 1].base[positionBuffer[j].offset + keyLength - k - 1] != (3^read.base[i + k]))
                        break;
                }
                positionBuffer[j].offset += i + keyLength - 1;
            }

            if (k == seedLength) {
                if (result[numResult].id == 0)
                    result[numResult] = positionBuffer[j];
                else
                    break;
            }
        }

        // check whther exist exact match kmer read or not
	  //seed not multi hit and not unmapped
        if (j == bufferSize && result[numResult].id != 0)
            ++numResult;
		//num of valid seed
    }
    if (numResult == 0) {
        result[0].id = 0;
        return result[0];
    }
    result.resize(numResult);
    sort(result.begin(), result.end());

    // select the highest overlapped kmer one
    i = firstMax = secondMax = result[numResult].id = 0;
    while (result[i].id != 0) {
        j = i;
        while (result[i].id == result[i+1].id && result[i].offset == result[i+1].offset)
            ++i;
        ++i;
        if (i - j >= firstMax) {
            secondMax = firstMax;
            firstMax = i - j;
            k = j;
        }
    }

    // if cannot decide top one, assign unmapped
    if (firstMax == secondMax) {
        result[0].id = 0;
        return result[0];
    }

    return result[k];
}



//////////////////////////////////////////////////////////////////////////////////////
// map each read and count the number of read mapped same contig and different one
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::mapPairMT(vector<SeqLib> &library, const long minInsertion, const long numThread)
{
    vector<platanus::Position> positionBuffer;
    platanus::SEQ forward;
    platanus::SEQ reverse;
    long totalPair = 0;
    long insertLength = 0;
    long numMappedSameContig = 0;
    long numMappedDiffContig = 0;
    long sum = 0;
    platanus::Position forwardPosition;
    platanus::Position reversePosition;
    char tmpChar;
    int i;

    cerr << "mapping reads..." << endl;

    omp_set_num_threads(numThread);


    #    pragma omp parallel for  schedule(static, 1) private(positionBuffer, forward, reverse, forwardPosition, reversePosition, insertLength) reduction(+: sum, totalPair, numMappedSameContig, numMappedDiffContig)
    for (i = 0; i < numThread; ++i) {
        if (library[i].insertLengthFP != NULL)
            fclose(library[i].insertLengthFP);
        library[i].insertLengthFP = platanus::makeTemporaryFile();

        if (library[i].mappedFP != NULL)
            fclose(library[i].mappedFP);
        library[i].mappedFP = platanus::makeTemporaryFile();
        rewind(library[i].pairFP);

        while (forward.readTemporaryFile(library[i].pairFP)) {
            reverse.readTemporaryFile(library[i].pairFP);
            ++totalPair;


            if (forward.length < seedLength || reverse.length < seedLength) continue;

            // find the contig where this read is mapped
            forwardPosition = this->mapRead(forward, positionBuffer);
            if (forwardPosition.id == 0) continue;
            reversePosition = this->mapRead(reverse, positionBuffer);
            if (reversePosition.id == 0) continue;
            // if forward read and reverse one are mapped same contig
            if (forwardPosition.id == -reversePosition.id) {
                if (forwardPosition.id > 0 && forwardPosition.offset < reversePosition.offset)
                    insertLength = static_cast<long>(reversePosition.offset - forwardPosition.offset + 1);
                else if (reversePosition.id > 0 && reversePosition.offset < forwardPosition.offset)
                    insertLength = static_cast<long>(forwardPosition.offset - reversePosition.offset + 1);
                else
                    continue;

                if (insertLength < std::min(forward.length, reverse.length)) continue;

                fwrite(&insertLength, sizeof(long), 1, library[i].insertLengthFP);
                ++numMappedSameContig;
            }
            else if (forwardPosition.id != reversePosition.id) {
                fwrite(&forwardPosition, sizeof(platanus::Position), 1, library[i].mappedFP);
                fwrite(&reversePosition, sizeof(platanus::Position), 1, library[i].mappedFP);
                ++numMappedDiffContig;
            }

            sum += forward.length + reverse.length;
        }
    }
    if (numMappedSameContig == 0) {
        throw platanus::MapError("No read mapped in the same contig!!");
    }
    if (numMappedDiffContig == 0) {
        throw platanus::MapError("no read links contigs!!");
    }

    library[0].setAverageCoverage(static_cast<double>(sum) / getSeqPoolSize());

    for (i = 1; i < numThread; ++i) {
        rewind(library[i].insertLengthFP);
        while (fread(&tmpChar, sizeof(char), 1, library[i].insertLengthFP))
            putc(tmpChar, library[0].insertLengthFP);
        fclose(library[i].insertLengthFP);
    }

    cerr << "TOTAL_PAIR = " << totalPair << endl;
    cerr << "MAPPED_PAIR = " << numMappedSameContig + numMappedDiffContig << " (" << static_cast<double>(numMappedSameContig + numMappedDiffContig) / totalPair << ")" << endl;
    cerr << "MAPPED_IN_DIFFERENT_CONTIGS = " << numMappedDiffContig << " (" << static_cast<double>(numMappedDiffContig) / totalPair << ")" << endl;
    cerr << "MAPPED_IN_SAME_CONTIG = " << numMappedSameContig << " (" << static_cast<double>(numMappedSameContig) / totalPair << ")" << endl;
    cerr << "AVERAGE_COVERAGE = " << library[0].getAverageCoverage() << endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// map single read which overlap gaps
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::mapSmallGap(vector<SeqLib> &library, FILE *gapFP, const long numThread)
{
    FILE *tmpFP[numThread];

    std::cerr << "mapping reads that cover small gaps..." << std::endl;

    omp_set_num_threads(numThread);

    # pragma omp parallel for schedule (static, 1) 
    for (long threadID = 0; threadID < numThread; ++threadID) {
        char gapSeq[platanus::ConstParam::MAX_READ_LEN];
        long gapSeqLength;
        vector<platanus::Position> leftBuffer;
        vector<platanus::Position> rightBuffer;
        platanus::Position leftResult;
        platanus::Position rightResult;
        platanus::SEQ read;
        Kmer31 leftKmer;
        Kmer31 rightKmer;

        tmpFP[threadID] = platanus::makeTemporaryFile();

        rewind(library[threadID].pairFP);
        while (read.readTemporaryFile(library[threadID].pairFP)) {
            if (read.length < 2 * this->seedLength) continue;
            long i;
            long start = 0;
            long end = 0;
            for (i = 0; i < this->keyLength; ++i) {
                if (read.base[i] == 4 || read.base[read.length - this->seedLength + i] == 4) break;
                leftKmer.forward = (leftKmer.forward << 2) | static_cast<unsigned long long>(read.base[i]);
                leftKmer.reverse = (leftKmer.reverse << 2) | static_cast<unsigned long long>((0x3^read.base[this->keyLength - i - 1]));
                rightKmer.forward = (rightKmer.forward << 2) | static_cast<unsigned long long>(read.base[read.length - this->seedLength + i]);
                rightKmer.reverse = (rightKmer.reverse << 2) | static_cast<unsigned long long>(0x3^read.base[read.length - this->seedLength + this->keyLength - i - 1]);
            }
            if (i != this->keyLength) continue;

            this->mapReadGapCloseLeft(read, leftKmer, leftBuffer);
            this->mapReadGapCloseRight(read, rightKmer, rightBuffer);


            leftResult.id = 0;
            rightResult.id = 0;

            for (i = 0; i < static_cast<long>(leftBuffer.size()); ++i) {
                unsigned long long j = 0;
                for (; j < rightBuffer.size(); ++j) {
                    if (leftBuffer[i].id == 0
                        || rightBuffer[j].id == 0
                        || leftBuffer[i].id != rightBuffer[j].id
                        || (leftBuffer[i].id > 0 && (leftBuffer[i].offset >= rightBuffer[j].offset || rightBuffer[j].offset - leftBuffer[i].offset + this->seedLength > 2 * read.length))
                        || (leftBuffer[i].id < 0 && (leftBuffer[i].offset <= rightBuffer[j].offset || leftBuffer[i].offset - rightBuffer[j].offset + this->seedLength > 2 * read.length))) continue;

                    if (leftResult.id != 0) break;
                    leftResult = leftBuffer[i];
                    rightResult = rightBuffer[j];
                }
                if (j != rightBuffer.size()) break;
            }
            if (static_cast<unsigned long long>(i) != leftBuffer.size() || leftResult.id == 0) continue;
            if (leftResult.id > 0) {
                start = searchGapStart(leftResult, rightResult, leftResult);
                if (start == 0) continue;
                end = searchGapEnd(leftResult, rightResult, leftResult, read.length);

                if (!checkGapContinuious(leftResult, rightResult, leftResult, start, end, read.length)) continue;

                gapSeqLength = 0;
                for (i = start; i < end; ++i) {
                    gapSeq[gapSeqLength] = read.base[i];
                    ++gapSeqLength;
                }
                leftResult.offset += start;
            } else {
                start = searchGapStart(rightResult, leftResult, leftResult);
                if (start == 0) continue;
                end = searchGapEnd(rightResult, leftResult, leftResult, read.length);

                if (!checkGapContinuious(rightResult, leftResult, leftResult, start, end, read.length)) continue;

                gapSeqLength = 0;
                for (i = read.length - start - 1; i >= read.length - end; --i) {
                    gapSeq[gapSeqLength] = 0x3 ^ read.base[i];
                    ++gapSeqLength;
                }
                leftResult.offset = rightResult.offset + start;
                leftResult.id *= -1;
            }

            gapSeqLength = end - start;
            fwrite(&leftResult, sizeof(platanus::Position), 1, tmpFP[threadID]);
            fwrite(&gapSeqLength, sizeof(long), 1, tmpFP[threadID]);
            if (gapSeqLength > 0)
                fwrite(gapSeq, sizeof(char), gapSeqLength, tmpFP[threadID]);
        }
    }

    char c;
    fseek(gapFP, 0, SEEK_END);
    for (long threadID = 0; threadID < numThread; ++threadID) {
        rewind(tmpFP[threadID]);
        while (fread(&c, sizeof(char), 1, tmpFP[threadID])) {
            putc(c, gapFP);
        }
        fclose(tmpFP[threadID]);
    }
}




//////////////////////////////////////////////////////////////////////////////////////
// wrapper function
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::mapReadGapCloseLeft(const platanus::SEQ &read, const Kmer31 &kmer, vector<platanus::Position> &buffer)
{
    this->mapReadGapClose(read, kmer, buffer, true);
}


//////////////////////////////////////////////////////////////////////////////////////
// wrapper function
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::mapReadGapCloseRight(const platanus::SEQ &read, const Kmer31 &kmer, vector<platanus::Position> &buffer)
{
    this->mapReadGapClose(read, kmer, buffer, false);
}


//////////////////////////////////////////////////////////////////////////////////////
// map single read actual for gap_close
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::mapReadGapClose(const platanus::SEQ &read, const Kmer31 &kmer, vector<platanus::Position> &buffer, const bool isLeft=true)
{
    long diff;
    long bufferSize = this->mapSeed(kmer, buffer);
    if (isLeft) {
        diff = 0;
    } else {
        diff = read.length - this->seedLength;
    }

    for (long i = 0; i < bufferSize; ++i) {
        long j = this->keyLength;
        if (buffer[i].id > 0) {
            if (buffer[i].offset > this->seq[buffer[i].id - 1].length - this->seedLength) continue;
            for (; j < this->seedLength; ++j) {
                if (this->seq[buffer[i].id - 1].base[buffer[i].offset + j] != read.base[diff + j]) break;
            }
        } else {
            if (buffer[i].offset < this->seedLength - this->keyLength) continue;
            for (; j < this->seedLength; ++j) {
                if (this->seq[-(buffer[i].id) - 1].base[buffer[i].offset + this->keyLength - j - 1] != (0x3^read.base[diff + j])) break;
            }
            buffer[i].offset -= this->seedLength - this->keyLength;
        }
        if (j != this->seedLength) {
            buffer[i].id = 0;
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////////////////////
long Mapper::searchGapStart(const platanus::Position &result1, const platanus::Position &result2, const platanus::Position &leftResult)
{
    long start = 0;
    for (long i = result1.offset + this->seedLength; i < result2.offset; ++i) {
        if (this->seq[abs(leftResult.id) - 1].base[i] == 4) {
            start = i - result1.offset;
            break;
        }
    }
    return start;
}



//////////////////////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////////////////////
long Mapper::searchGapEnd(const platanus::Position &result1, const platanus::Position &result2, const platanus::Position &leftResult, const long length)
{
    long end = 0;
    for (long i = result2.offset - 1; i >= result1.offset + this->seedLength; --i) {
        if (this->seq[abs(leftResult.id) - 1].base[i] == 4) {
            end = length - (result2.offset + this->seedLength - 1 - i);
            break;
        }
    }
    return end;
}



//////////////////////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////////////////////
bool Mapper::checkGapContinuious(const platanus::Position &result1, const platanus::Position &result2, const platanus::Position &leftResult, const long start, const long end, const long length)
{
    long i;
    for (i = result1.offset + start; i < result2.offset - (length - end - this->seedLength); ++i) {
        if (this->seq[abs(leftResult.id) - 1].base[i] != 4) break;
    }
    return i == result2.offset - (length - end - this->seedLength);
}



//////////////////////////////////////////////////////////////////////////////////////
// map each read and count the number of read mapped same contig
// for gap close function
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::gatherPairReadMappedSameContig(vector<SeqLib> &library, const long numThread)
{
    vector<platanus::Position> positionBuffer;
    platanus::SEQ forward;
    platanus::SEQ reverse;
    long totalPair = 0;
    long insertLength = 0;
    long numMappedSameContig = 0;
    platanus::Position forwardPosition;
    platanus::Position reversePosition;
    char tmpChar;
    int i;

    cerr << "mapping reads..." << endl;

    omp_set_num_threads(numThread);


    #    pragma omp parallel for  schedule(static, 1) private(positionBuffer, forward, reverse, forwardPosition, reversePosition, insertLength) reduction(+: totalPair, numMappedSameContig)
    for (i = 0; i < numThread; ++i) {
        if (library[i].insertLengthFP != NULL)
            fclose(library[i].insertLengthFP);
        library[i].insertLengthFP = platanus::makeTemporaryFile();

        if (library[i].mappedFP != NULL)
            fclose(library[i].mappedFP);
        library[i].mappedFP = platanus::makeTemporaryFile();
        rewind(library[i].pairFP);

        while (forward.readTemporaryFile(library[i].pairFP)) {
            reverse.readTemporaryFile(library[i].pairFP);
            ++totalPair;

            if (forward.length < seedLength || reverse.length < seedLength) continue;

            forwardPosition = this->mapRead(forward, positionBuffer);
            reversePosition = this->mapRead(reverse, positionBuffer);
            if (forwardPosition.id == 0 && reversePosition.id == 0) continue;

            fwrite(&forwardPosition, sizeof(platanus::Position), 1, library[i].mappedFP);
            for (int j = 0; j < forward.numUnknown; ++j)
                forward.base[forward.positionUnknown[j]] = 4;
            fwrite(&(forward.length), sizeof(long), 1, library[i].mappedFP);
            fwrite(forward.base.c_str(), sizeof(char), forward.length, library[i].mappedFP);

            fwrite(&reversePosition, sizeof(platanus::Position), 1, library[i].mappedFP);
            for (int j = 0; j < reverse.numUnknown; ++j)
                reverse.base[reverse.positionUnknown[j]] = 4;
            fwrite(&(reverse.length), sizeof(long), 1, library[i].mappedFP);
            fwrite(reverse.base.c_str(), sizeof(char), reverse.length, library[i].mappedFP);
            
            if (forwardPosition.id == -reversePosition.id) {
                if (forwardPosition.id > 0 && forwardPosition.offset < reversePosition.offset) {
                    insertLength = static_cast<long>(reversePosition.offset) - forwardPosition.offset;
                } else if (reversePosition.id > 0 && forwardPosition.offset > reversePosition.offset) {
                    insertLength = static_cast<long>(forwardPosition.offset) - reversePosition.offset;
                } else
                    continue;
                
                if (insertLength < std::min(forward.length, reverse.length)) continue;
                fwrite(&insertLength, sizeof(long), 1, library[i].insertLengthFP);
                ++numMappedSameContig;
            }
        }
    }

    if (numMappedSameContig == 0) {
        throw platanus::MapError("no read mapped in the same contig!!");
    }

    for (i = 1; i < numThread; ++i) {
        rewind(library[i].insertLengthFP);
        while (fread(&tmpChar, sizeof(char), 1, library[i].insertLengthFP))
            putc(tmpChar, library[0].insertLengthFP);
        fclose(library[i].insertLengthFP);
    }

    cerr << "TOTAL_PAIR = " << totalPair << endl;
    cerr << "MAPPED_IN_SAME_CONTIG = " << numMappedSameContig << " (" << static_cast<double>(numMappedSameContig) / totalPair << ")" << endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// decide the contig where the variable read is mapped
// contig contains bubble seq
//////////////////////////////////////////////////////////////////////////////////////
platanus::Position HeteroMapper::mapRead(const platanus::SEQ &read, vector<platanus::Position> &positionBuffer)
{
    long numResult = 0;
    vector<platanus::Position> result(read.length/seedLength + 2);
    Kmer31 kmer;
    long firstMax, secondMax;
    int i, j, k = 0;
    for (i = read.length - seedLength; i > -(seedLength); i -= seedLength) {
        if (i < 0)
            i = 0;
        kmer.forward = kmer.reverse = 0;
        for (j = 0; j < keyLength; ++j) {
            if (read.base[i + j] == 4) break;
            kmer.forward = (kmer.forward << 2) | static_cast<unsigned long long>(read.base[i + j]);
            kmer.reverse = (kmer.reverse << 2) | static_cast<unsigned long long>(0x3^read.base[i + keyLength - j - 1]);
        }
        if (j != keyLength) continue;
        // make this function
        contigMap.mapKey(kmer, positionBuffer);
        result[numResult].id = 0;
        auto it = positionBuffer.begin();
        auto end = positionBuffer.end();
        for (; it != end; ++it) {
            if (it->id > 0) {
                if (it->offset > contigMap.getSeqLength(it->id - 1) - seedLength)
                    continue;
                for (k = keyLength; k < seedLength; ++k) {
                    if (contigMap.getBase(it->id - 1, it->offset + k) != read.base[i + k])
                        break;
                }
                it->offset -= i;
            }
            else {
                if (it->offset < seedLength - keyLength)
                    continue;
                for (k = keyLength; k < seedLength; ++k) {
                    if (contigMap.getBase(-(it->id) - 1, it->offset + keyLength - k - 1) != (3^read.base[i + k]))
                        break;
                }
                it->offset += i + keyLength - 1;
            }

            if (k == seedLength) {
                if (result[numResult].id == 0)
                    result[numResult] = *it;
                else
                    break;
            }
        }
        if (it == end) {
            if (result[numResult].id != 0) {
                ++numResult;
                continue;
            }
        }
        else
            continue;

        bubbleMap.mapKey(kmer, positionBuffer);
        result[numResult].id = 0;

        it = positionBuffer.begin();
        end = positionBuffer.end();

        for (; it != end; ++it) {
            if (it->id > 0) {
                if (it->offset > bubbleMap.getSeqLength(it->id - 1) - seedLength)
                    continue;
                for (k = keyLength; k < seedLength; ++k) {
                    if (bubbleMap.getBase(it->id - 1, it->offset + k) != read.base[i + k])
                        break;
                }
                it->offset -= i;
                if (bubblePosition[it->id - 1].id > 0) {
                    it->offset = bubblePosition[it->id - 1].offset + it->offset - bubbleMap.getSeqLength(it->id - 1) / 2;
                    it->id = bubblePosition[it->id - 1].id;
                }
                else {
                    it->offset = bubblePosition[it->id - 1].offset - it->offset + bubbleMap.getSeqLength(it->id - 1) / 2;
                    it->id = bubblePosition[it->id - 1].id;
                }
            }
            else {
                if (it->offset < seedLength - keyLength)
                    continue;
                for (k = keyLength; k < seedLength; ++k) {
                    if (bubbleMap.getBase(-(it->id) - 1, it->offset + keyLength - k - 1) != (0x3^read.base[i + k]))
                        break;
                }
                it->offset += i + keyLength - 1;
                if (bubblePosition[-(it->id) - 1].id > 0) {
                    it->offset = bubblePosition[-(it->id) - 1].offset + it->offset - bubbleMap.getSeqLength(-(it->id) - 1) / 2;
                    it->id = -(bubblePosition[-(it->id) - 1].id);
                }
                else {
                    it->offset = bubblePosition[-(it->id) - 1].offset - it->offset + bubbleMap.getSeqLength(-(it->id) - 1) / 2;
                    it->id = -(bubblePosition[-(it->id) - 1].id);
                }
            }

            if (k == seedLength) {
                if (result[numResult].id == 0)
                    result[numResult] = *it;
                else
                    break;
            }
        }
        if (it == end && result[numResult].id != 0)
            ++numResult;
    }
    if (numResult == 0) {
        result[0].id = 0;
        return result[0];
    }
    result.resize(numResult + 1);
    sort(result.begin(), result.end() - 1);

    i = firstMax = secondMax = result[numResult].id = 0;
    while (result[i].id != 0) {
        j = i;
        while (result[i].id == result[i+1].id && result[i].offset == result[i+1].offset)
            ++i;
        ++i;
        if (i - j >= firstMax) {
            secondMax = firstMax;
            firstMax = i - j;
            k = j;
        }
    }

    if (firstMax == secondMax) {
        result[0].id = 0;
        return result[0];
    }

    return result[k];
}


//////////////////////////////////////////////////////////////////////////////////////
// map each read and count the number of read mapped same contig (contain bubbles) and different one
//////////////////////////////////////////////////////////////////////////////////////
void HeteroMapper::mapPairMT(vector<SeqLib> &library, const long minInsertion, const long numThread)
{
    vector<platanus::Position> positionBuffer;
    platanus::SEQ forward;
    platanus::SEQ reverse;
    long totalPair = 0;
    long insertLength = 0;
    long numMappedSameContig = 0;
    long numMappedDiffContig = 0;
    long sum = 0;
    platanus::Position forwardPosition;
    platanus::Position reversePosition;
    char tmpChar;
    int i;
    if (bubbleMap.getNumSeq() == 0) {
        contigMap.mapPairMT(library, minInsertion, numThread);
        return;
    }

    cerr << "mapping reads..." << endl;

    omp_set_num_threads(numThread);
    #    pragma omp parallel for  schedule(static, 1) private(positionBuffer, forward, reverse, forwardPosition, reversePosition, insertLength) reduction(+: sum, totalPair, numMappedSameContig, numMappedDiffContig)
    for (i = 0; i < numThread; ++i) {

        if (library[i].insertLengthFP != NULL)
            fclose(library[i].insertLengthFP);
        library[i].insertLengthFP = platanus::makeTemporaryFile();

        if (library[i].mappedFP != NULL)
            fclose(library[i].mappedFP);
        library[i].mappedFP = platanus::makeTemporaryFile();
        rewind(library[i].pairFP);

        while (forward.readTemporaryFile(library[i].pairFP)) {
            reverse.readTemporaryFile(library[i].pairFP);
            ++totalPair;
            if (forward.length < seedLength || reverse.length < seedLength) continue;
            forwardPosition = contigMap.mapRead(forward, positionBuffer);
            if (forwardPosition.id == 0) {
                forwardPosition = this->mapRead(forward, positionBuffer);
                if (forwardPosition.id == 0) continue;
            }

            reversePosition = contigMap.mapRead(reverse, positionBuffer);
            if (reversePosition.id == 0) {
                reversePosition = this->mapRead(reverse, positionBuffer);
                if (reversePosition.id == 0) continue;
            }


            if (forwardPosition.id == 0 || reversePosition.id == 0) continue;
            if (forwardPosition.id == -reversePosition.id) {
                if (forwardPosition.id > 0 && forwardPosition.offset < reversePosition.offset)
                    insertLength = static_cast<long>(reversePosition.offset - forwardPosition.offset + 1);
                else if (reversePosition.id > 0 && reversePosition.offset < forwardPosition.offset)
                    insertLength = static_cast<long>(forwardPosition.offset - reversePosition.offset + 1);
                else
                    continue;

                if (insertLength < minInsertion || insertLength < std::min(forward.length, reverse.length))
                    continue;

                fwrite(&insertLength, sizeof(long), 1, library[i].insertLengthFP);
                ++numMappedSameContig;
            }
            else if (forwardPosition.id != reversePosition.id) {
                fwrite(&forwardPosition, sizeof(platanus::Position), 1, library[i].mappedFP);
                fwrite(&reversePosition, sizeof(platanus::Position), 1, library[i].mappedFP);
                ++numMappedDiffContig;
            }
            sum += forward.length + reverse.length;
        }
    }

    if (numMappedSameContig == 0) {
        throw platanus::MapError("no read mapped in the same contig!!");
    }
    if (numMappedDiffContig == 0) {
        throw platanus::MapError("no read links contigs!!");
    }

    library[0].setAverageCoverage(static_cast<double>(sum) / contigMap.getSeqPoolSize());
    for (i = 1; i < numThread; ++i) {
        rewind(library[i].insertLengthFP);
        while (fread(&tmpChar, sizeof(char), 1, library[i].insertLengthFP))
            putc(tmpChar, library[0].insertLengthFP);
        fclose(library[i].insertLengthFP);
    }
    cerr << "TOTAL_PAIR = " << totalPair << endl;
    cerr << "MAPPED_PAIR = " << numMappedSameContig + numMappedDiffContig << " (" << static_cast<double>(numMappedSameContig + numMappedDiffContig) / totalPair << ")" << endl;
    cerr << "MAPPED_IN_DIFFERENT_CONTIGS = " << numMappedDiffContig << " (" << static_cast<double>(numMappedDiffContig) / totalPair << ")" << endl;
    cerr << "MAPPED_IN_SAME_CONTIG = " << numMappedSameContig << " (" << static_cast<double>(numMappedSameContig) / totalPair << ")" << endl;
    cerr << "AVERAGE_COVERAGE = " << library[0].getAverageCoverage() << endl;
}

