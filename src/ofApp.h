#pragma once

#include "ofMain.h"
#include "ofxFft.h"
#include "ofxOsc.h"

// listening port
#define PORT 8800
#define CH 2
// max number of strings to display
#define NUM_MSG_STRINGS 20

enum {SINE, MIC, NOISE};

typedef struct {
    
    float     x;
    float     y;
    bool     bBeingDragged;
    bool     bOver;
    float     radius;
    
}draggableVertex;

class ofApp : public ofBaseApp {
public:
    
    // OF //
	void setup();
    void update();
    void draw();
    
    // FFT //
	void plot(vector<float>& buffer, float scale, float offset);
	void audioReceived(float* input, int bufferSize, int nChannels);
	int plotHeight, bufferSize;
	int spectrogramOffset;
	int mode;
	int appWidth, appHeight;
    ofImage spectrogram;
    ofxFft* fft;
	ofMutex soundMutex;
	vector<float> drawBins, middleBins, audioBins;
	vector<float> drawBuffer, middleBuffer, audioBuffer;
    
    // Events //
    void keyPressed(int key);
    void keyReleased(int key);
    void mouseMoved(int x, int y);
    void mouseDragged(int x, int y, int button);
    void mousePressed(int x, int y, int button);
    void mouseReleased(int x, int y, int button);
    void mouseEntered(int x, int y);
    void mouseExited(int x, int y);
    void windowResized(int w, int h);
    void dragEvent(ofDragInfo dragInfo);
    void gotMessage(ofMessage msg);
    
    // OSC //
    int mouseButtonInt = 0;
    int currentMsgString;
    string msgStrings[NUM_MSG_STRINGS];
    string mouseButtonState = "";
    float timers[NUM_MSG_STRINGS];
    float mouseXf = 0;
    float mouseYf = 0;
    ofTrueTypeFont font;
    ofxOscReceiver receiver;
    ofImage receivedImage;
    
    // VIDEO //
    ofVideoPlayer vid;
    bool frameByframe = true;
    
    // Graphics //
    ofVec3f previous, current;
    ofEasyCam easyCam;
    
	// Shader //
	ofShader shader;

    deque<ofVec3f> pathVertices;
    ofMesh pathLines;
    
    
};
