#define _USE_MATH_DEFINES

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;


vector<vector<double>> window(string data, int frameLength, int hopSize, string windowType){

	vector<vector<double>>	frameMatrix = {};
	vector<double> 			frame = {};
	vector<double> 			windowFunction = {};	
	vector<double> 			dataV = {};

	ifstream inFile;
	inFile.open("data/list.dat");
	ofstream oFile;
	oFile.open("data/wlist.dat");

	double x;
	while(inFile >> x){
		dataV.push_back(x);
	}													//All data moved from list.dat to dataV, kinda unnecessary to put it in textfile, could fix
	inFile.close();
	int numberOfFrames = 1 + floor((dataV.size() - frameLength) / hopSize);
		
	if(windowType == "hamming"){						//Implement Hamming window  w(n)=0.54−0.46cos(2πn/M−1),   0≤n≤M−1
		for(int i = 0; i < frameLength - 1; i++){ 
			double point = 0.54 - 0.46*cos((2*M_PI*i)/(frameLength-1));
			windowFunction.push_back(point);
		}						
	}
	else cout << "Windowing function not supported" << endl;

	for(int i = 0; i < numberOfFrames; i++){ 						//for each row
		for(int j = 0; j < frameLength; j++){						//for each column
			double dataPoint = dataV[i*hopSize + j];
			if(dataPoint){
				dataPoint *= windowFunction[j];						//assign correct data point * windowing function point
				frame.push_back(dataPoint); 						
				oFile << dataPoint << endl;
			} else oFile << 0 << endl;
		}
		frameMatrix.push_back(frame);
		frame.clear();
	}
	oFile.close();
	return frameMatrix;
}



