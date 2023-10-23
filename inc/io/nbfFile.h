#ifndef FILE_nbfFile
#define FILE_nbfFile

class nbfFile
{

public:

	nbfFile(){
		this->fileName = NULL;
	}

	// Set file name.
	void setFileName( const char * _fileName ){
		if (this->fileName){
			delete [] this->fileName;
		}
		if (_fileName){
			this->fileName = new char[strlen(_fileName)+1];
			strcpy( this->fileName,_fileName );
		}
		else{
			this->fileName = NULL;
		}
	}

protected:

	// File name.
	char * fileName;

};


#endif // FILE_nbfFile