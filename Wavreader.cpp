#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string.h>

using namespace std;
									//double byteToDouble( char firstByte, char secondByte );
									//WAVE PCM soundfile format (you can find more in https://ccrma.stanford.edu/courses/422/projects/WaveFormat/)

struct wavHeader {			
    unsigned char 	chunk_id[4];
    int 			chunk_size;
    unsigned char 	audioFormat[4];
    unsigned char 	subchunk1_id[4];
    int 			subchunk1_size;
    unsigned char 	audio_format[2];
    unsigned char 	num_channels[2];
    int 			sample_rate;			//sample_rate denotes the sampling rate.
    int 			byte_rate;
    unsigned char	block_align[2];
    unsigned char	bits_per_sample[2];
    unsigned char 	subchunk2_id[4];
    int				subchunk2_size; 		//subchunk2_size denotes the number of samples.
    										//char data; // actual data : Added by tarmizi
};

int readWav(){

    wavHeader* newWav = (wavHeader*)malloc(sizeof(wavHeader) + 1);

	FILE* myFile;							//Open wave file in read mode
	myFile = fopen("Sounds/SX83.wav","rb");
	ofstream oFile;
	oFile.open("data/list.dat");

	fseek (myFile, 0, SEEK_END);   			// non-portable
    int fileSize = ftell(myFile); 			//Returns the current value of the position indicator of the stream.
                                			//For binary streams, this is the number of bytes from the beginning of the file.
    rewind(myFile); 						//Sets the position indicator associated with stream to the beginning of the file.

    cout << "\t**********Header information********" << endl;			//print file information
    cout << "The size of the file is " << fileSize << " bytes." << endl << endl;

	fread(newWav->chunk_id,1,4,myFile);
	cout << newWav->chunk_id << endl;

	fread(&newWav->chunk_size,1,4,myFile);
	cout << newWav->chunk_size << endl;

	fread(newWav->audio_format,1,4,myFile);
    cout << newWav->audio_format << endl;

    fread(newWav->subchunk1_id,1,4,myFile);
    cout << newWav->subchunk1_id << endl;

    fread(&newWav->subchunk1_size,1,4,myFile);
    cout << newWav->subchunk1_size << endl;

    fread(newWav->audioFormat,1,2,myFile);
    cout << (int)newWav->audioFormat[0] << endl;

    fread(newWav->num_channels,1,2,myFile);
    int num_channels = (int)*newWav->num_channels;
    cout <<"Number of channels: "<< num_channels << endl;

    fread(&newWav->sample_rate,1,4,myFile);
    cout <<"Sampling rate is: "<< newWav->sample_rate << endl;

    fread(&newWav->byte_rate,1,4,myFile);
    cout <<"Byte rate is: "<< newWav->byte_rate << endl;

    fread(newWav->block_align,1,2,myFile);
    cout << (int)*newWav->block_align << endl;

    fread(newWav->bits_per_sample,1,2,myFile);
    int bits_per_sample = (int)*newWav->bits_per_sample;
    cout << bits_per_sample << endl;

    if(newWav->subchunk1_size > 16){
        char dummy[newWav->subchunk1_size - 16];
        fread(dummy, 1, newWav->subchunk1_size-16, myFile);
    }

    fread(newWav->subchunk2_id, 1, 4, myFile);
    cout << newWav->subchunk2_id << endl;

    fread(&newWav->subchunk2_size, 1, 4, myFile);
    cout << newWav->subchunk2_size << endl;

    if(bits_per_sample == 16){
    	cout << "Number of Samples " << newWav->subchunk2_size/(num_channels*(bits_per_sample/8)) << endl << endl;
	}

    cout << "The size of the wav header is " << (sizeof(wavHeader) + newWav->subchunk1_size - 16) << " bytes" << endl;
    if(bits_per_sample == 16){
   		int numOfSamples = newWav->subchunk2_size/(num_channels*(bits_per_sample/8));
       // int *data = (int*)malloc(sizeof(int)*numOfSamples);
        int counter = 0;
        short int readBit0;
        short int readBit1;
        //double conv0;
        //double conv1;

		while(counter <= numOfSamples){
            if(num_channels == 2){
            	fread(&readBit0, 1, 2, myFile);
                fread(&readBit1, 1, 2, myFile);
            }
            if(num_channels == 1) fread(&readBit1, 1, 2, myFile);

            //conv0 = readBit0 / 32768.0; // Converting the sample values between -1 and 1 for plotting a waveform
            //conv1 = readBit1 / 32768.0;  						
           // oFile << readBit0 << endl;//conv0 << endl; // Write the sample values to a txt file.
            oFile << readBit1 << endl;//conv1 << endl; //
            counter++;
		}
		fclose(myFile);
        cout << "Number of samples actually read in " << counter - 1 << endl;
	}
	int sampleRate = newWav->sample_rate;
	oFile.close();
    return sampleRate;
}

