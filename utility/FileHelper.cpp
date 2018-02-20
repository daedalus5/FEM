
#include <fstream>
#include <iostream>
#include "FileHelper.h"

bool FileHelper::readFloats(char *path, std::vector<float> &result) {

    std::ifstream inFile;
    float value;

    inFile.open(path);

    if (!inFile)
    {
        std::cout << "\nError opening file.\n";
        return false;
    }

    while (inFile >> value)
    {
        result.push_back(value);
    }

    inFile.close();

    return true;
}

bool FileHelper::printFloats(char *path, std::vector<float> &output) {


    std::ofstream outputStream;

    outputStream.open(path);

    if (!outputStream)
    {
        std::cout << "\nError writing to file.\n";
        return false;
    }

    for(int i = 0; i < output.size(); ++i)
    {
        outputStream << output[i] << " ";
    }

    outputStream.close();

    return true;
}