#include "matplotlibcpp.h"
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "string.h"
#include "windowing.cpp"
#include "Wavreader.cpp"
#include <fftw3.h>
#include "resample.h"

using namespace std;

int main(){

	//1.1 Read the audio file and sampling rate
	int sampleRate = readWav();
	int sampleRateTarget = 16000;
	cout << "Sample rate is: " << sampleRate << endl;

	vector<double> dataV = {};
	ifstream inFile;
	inFile.open("data/list.dat");
	double tp;
	while(inFile >> tp){
		dataV.push_back(tp);
	}											//All data moved from list.dat to dataV, kinda unnecessary to put it in textfile, could fix
	inFile.close();

	//1.2 Make sure the sampling rate is 16kHz, resample if not
	vector<double> dataResampled = {};
	if(sampleRate != sampleRateTarget){
		resample(sampleRateTarget,sampleRate,dataV,dataResampled);
	}

	//1.3 Split the data sequence into windows, implement windowing.cpp
	int frameLength = sampleRate * 0.025;  //25ms in samples, (400?)
	int hopSize 	= sampleRate * 0.0125; //12.5 ms in samples (50% overlap)
	//vector<string> windowTypes = ("rect", "hann", "cosine", "hamming");
	vector<vector<double>> frameMatrix = window("data/list.dat", frameLength, hopSize, "hamming"); //Windowing

	

	//TODO 1.4 Visualization?
	namespace plt = matplotlibcpp;
	using namespace plt;

    figure(1);

	subplot(3,1,1);
	vector<double> t_axis1 = {};  //divide 0 - length of dataV by 16000
	int size = dataV.size();
	for(int i = 0; i < size; i++){
		t_axis1.push_back(i/16000);
	}
	plot(t_axis1, dataV);
 	xlim(0,2);


	subplot(3,1,2);
	vector<double> t_axis2 = {}; //divide 0 - frameLength with 16
	for(float i = 0; i < frameLength; i++){
		t_axis2.push_back(i/(sampleRate/1000));
	}
	vector<double> temp = {};
	unsigned int count1 = 0;
	for(int i = 0; i < frameLength; i++){
		temp.push_back(frameMatrix[22][i]); //frameMatrix[:,23]?
		count1++;
		if(count1 == t_axis2.size()) break;
	}
	plot(t_axis2, temp); 
	ylim(-1500, 2000);

	//FFTW
	double matrixSize = frameMatrix.size()*frameLength;					//total size of frameMatrix (n x m)
	double* in = (double*) fftw_malloc(sizeof(double) * matrixSize);	//allocate fft-input array
	fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * matrixSize/2 + 1);	//allocate fft-output array
    fftw_plan p;	//define fft plan
	int N = 1024;
	ofstream offtfile;
	offtfile.open("data/fftdata.txt"); 	//open textfile for data to matlab for testing
	vector<double> outV = {};		//vector for fft-output data
    p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE); //init fft-plan
	for(unsigned int i = 0; i < frameMatrix.size(); i++){		//for each column in frameMatrix, set it as fft-input
		for(int j = 0; j < frameLength; j++){			//
			in[j] = frameMatrix[i][j];					//
		}												//
		fftw_execute(p);								//execute fft for each column
		for(int i = N/2+1; i > 0; i--){												//
			offtfile << 20*log10(abs(out[i][0] + out[i][1])) << std::endl;								//send output to textfile
			outV.push_back(20*log10(abs(out[i][0] + out[i][1])));	//push output to outV - datavector
		}
	}    
	offtfile.close();		//close outputfile
	fftw_destroy_plan(p);	//destruct fft_plan
    fftw_free(in); 			//free allocated memory
	fftw_free(out);			//

	ofstream oenerfile;
	oenerfile.open("data/energydata.txt");
	unsigned int nEnerBands = 32;
	vector<vector<double>> enerOut = {};
	vector<double> enerFrame = {};
	//apply ener_bands to outV
	int nfft = 216;//outV.column.size(), 216
	//int nframes = 513; //outV.row.size()
	int bandsPerEner = floor(nfft/nEnerBands);
	cout << bandsPerEner << endl;
	unsigned int count = 0;
	int sum = 0;
	int mean = 0;
	
	for(auto i = outV.begin(); i != outV.end(); i++){
		sum += pow(abs(*i),2);
		
		if(count % bandsPerEner == 0 /*&& count != 0*/){
			enerFrame.push_back(10*log10(sum));
			oenerfile << 10*log10(sum) << endl;
			mean += 10*log10(sum);
			sum = 0;
		//	cout << "joo " << *i << endl;
		}				
		if(enerFrame.size() == nEnerBands){
			mean /= enerFrame.size();
			for(unsigned int k = 0; k < enerFrame.size(); k++){			
				enerFrame[k] -= mean;		//removing mean?
			}
			enerOut.push_back(enerFrame);
			enerFrame.clear();
		}		
		count++;	
	}
	oenerfile.close();

	subplot(3,1,3);
	vector<double> f_axis = {};
	for(unsigned int i = 0; i < (N/2+1)*frameMatrix.size(); i++){
		f_axis.push_back(sampleRate * i / N);		
	}
	plot(f_axis,outV);
	xlim(-0.05,0.05);
	ylim(0,100);
    show();

	vector<vector<double>> frameMatrixFft = {};			
	vector<double> fmfFrame = {};
	for(unsigned int i = 0; i < frameMatrix.size(); i++){
		for(int j = 0; j < N/2+1; j++){
			fmfFrame.push_back(outV[j + i*(N/2+1)]);
		}
		frameMatrixFft.push_back(fmfFrame);
		fmfFrame.clear();
	}
	vector<vector<double>> x,y,z = {};
	
	for(unsigned int i = 0; i < frameMatrix.size(); i+=10){
		vector<double> xrow, yrow, zrow = {};
		for(int j = 0; j < N/2+1; j+=10){
			xrow.push_back(i);
            yrow.push_back(j);
            zrow.push_back(frameMatrixFft[i][j]);
		}
		x.push_back(xrow);
        y.push_back(yrow);
        z.push_back(zrow);
	}
	/*for (double i = -5; i <= 5;  i += 0.25) {
        std::vector<double> xrow, yrow, zrow;
        for (double j = -5; j <= 5; j += 0.25) {
            xrow.push_back(i);
            yrow.push_back(j);
            zrow.push_back(::std::sin(::std::hypot(i, j)));
        }
        x.push_back(xrow);
        y.push_back(yrow);
        z.push_back(zrow);
    }*/

	plot_surface(x,y,z);
	show();


	return 0;
}
