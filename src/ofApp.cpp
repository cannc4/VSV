#include "ofApp.h"

float powFreq(float i) {
    return powf(i, 3);
}

void rotateToNormal(ofVec3f normal) {
    normal.normalize();
    
    float rotationAmount;
    ofVec3f rotationAngle;
    ofQuaternion rotation;
    
    ofVec3f axis(0, 0, 1);
    rotation.makeRotate(axis, normal);
    rotation.getRotate(rotationAmount, rotationAngle);
    ofRotateDeg(rotationAmount, rotationAngle.x, rotationAngle.y, rotationAngle.z);
}

void ofApp::setup() {
    // Graphics init //
    ofSetVerticalSync(true);
    // draw the vertices in pathLines as a line strip
    pathLines.setMode(OF_PRIMITIVE_LINE_STRIP);
    
    // Video init //
    vid.load("vid/test.mov");
    vid.setLoopState(OF_LOOP_NORMAL);
    vid.play();
    ofSetVerticalSync(true);
    frameByframe = false;

    // OSC init //
    // listen on the given port
    ofLog() << "listening for osc messages on port " << PORT;
    receiver.setup(PORT);
	
    // FFT init //
	plotHeight = 128;
	bufferSize = 512;

	fft = ofxFft::create(bufferSize, OF_FFT_WINDOW_HAMMING);
	// To use FFTW, try:
	//fft = ofxFft::create(bufferSize, OF_FFT_WINDOW_HAMMING, OF_FFT_FFTW);

	spectrogram.allocate(bufferSize, fft->getBinSize(), OF_IMAGE_GRAYSCALE);
    spectrogram.setColor(ofColor::black);
    spectrogramOffset = 0;

	drawBuffer.resize(bufferSize);
	middleBuffer.resize(bufferSize);
	audioBuffer.resize(bufferSize);
	
	drawBins.resize(fft->getBinSize());
	middleBins.resize(fft->getBinSize());
	audioBins.resize(fft->getBinSize());

	// 0 output channels,
	// 1 input channel
	// 44100 samples per second
	// [bins] samples per buffer
	// 4 num buffers (latency)

	ofSoundStreamSetup(0, CH, this, 44100, bufferSize, 4);

	mode = SINE;
	appWidth = ofGetWidth();
	appHeight = ofGetHeight();

	ofBackground(0, 0, 0);
}

void ofApp::update(){
    vid.update();
    // hide old messages
    for(int i = 0; i < NUM_MSG_STRINGS; i++){
        if(timers[i] < ofGetElapsedTimef()){
            msgStrings[i] = "";
        }
    }
    
    // check for waiting messages
    while(receiver.hasWaitingMessages()){
        
        // get the next message
        ofxOscMessage m;
        receiver.getNextMessage(m);
        
        // check for mouse moved message
        if(m.getAddress() == "/mouse/position"){
            
            // both the arguments are floats
            mouseXf = m.getArgAsFloat(0);
            mouseYf = m.getArgAsFloat(1);
        }
        // check for mouse button message
        else if(m.getAddress() == "/mouse/button"){
            
            // first argument is int32, second is a string
            mouseButtonInt = m.getArgAsInt32(0);
            mouseButtonState = m.getArgAsString(1);
        }
        // check for an image being sent
        // note: the size of the image depends greatly on your network buffer
        // sizes, if an image is too big the message won't come through
        else if(m.getAddress() == "/image"){
            ofBuffer buffer = m.getArgAsBlob(0);
            receivedImage.load(buffer);
        }
        else{
            
            // unrecognized message: display on the bottom of the screen
            string msgString;
            msgString = m.getAddress();
            msgString += ":";
            for(size_t i = 0; i < m.getNumArgs(); i++){
                
                // get the argument type
                msgString += " ";
                msgString += m.getArgTypeName(i);
                msgString += ":";
                
                // display the argument - make sure we get the right type
                if(m.getArgType(i) == OFXOSC_TYPE_INT32){
                    msgString += ofToString(m.getArgAsInt32(i));
                }
                else if(m.getArgType(i) == OFXOSC_TYPE_FLOAT){
                    msgString += ofToString(m.getArgAsFloat(i));
                }
                else if(m.getArgType(i) == OFXOSC_TYPE_STRING){
                    msgString += m.getArgAsString(i);
                }
                else{
                    msgString += "unhandled argument type " + m.getArgTypeName(i);
                }
            }
            
            // add to the list of strings to display
            msgStrings[currentMsgString] = msgString;
            timers[currentMsgString] = ofGetElapsedTimef() + 5.0f;
            currentMsgString = (currentMsgString + 1) % NUM_MSG_STRINGS;
            
            // clear the next line
            msgStrings[currentMsgString] = "";
        }
    }
    
    // Graphics Update //
    previous = current;
    
    // generate a noisy 3d position over time
    float t = (2 + ofGetElapsedTimef()) * .1;
    current.x = ofSignedNoise(t, 0, 0);
    current.y = ofSignedNoise(0, t, 0);
    current.z = ofSignedNoise(0, 0, t);
    current *= 400; // scale from -1,+1 range to -400,+400
    
    // add the current position to the pathVertices deque
    pathVertices.push_back(current);
    // if we have too many vertices in the deque, get rid of the oldest ones
    while(pathVertices.size() > 200) {
        pathVertices.pop_front();
    }
    
    // clear the pathLines ofMesh from any old vertices
    pathLines.clear();
    // add all the vertices from pathVertices
    for(unsigned int i = 0; i < pathVertices.size(); i++) {
        pathLines.addVertex(pathVertices[i]);
    }
}

void ofApp::draw() {
    
    // Graphics //
    ofColor cyan = ofColor::fromHex(0x00abec);
    ofColor magenta = ofColor::fromHex(0xec008c);
    ofColor yellow = ofColor::fromHex(0xffee00);
    
    ofBackgroundGradient(magenta * .6, magenta * .4);
    ofNoFill();
    
    easyCam.begin();
    ofRotateXDeg(15);
    
    ofSetColor(0);
    ofDrawGrid(500, 10, false, false, true, false);
    
    // draw the path of the box
    ofSetLineWidth(2);
    ofSetColor(cyan);
    pathLines.draw();
    
    // draw a line connecting the box to the grid
    ofSetColor(yellow);
    ofDrawLine(current.x, current.y, current.z, current.x, 0, current.z);
    
    ofTranslate(current.x, current.y, current.z);
    if( (current - previous ).length() > 0.0 ){
        // translate and rotate every 3D object after this line to the current position and orientation of our line, but only if the line is longer than 0 or has a length
        rotateToNormal(current - previous);
    }
    ofSetColor(255);
    ofDrawBox(32);
    ofDrawAxis(32);
    
    easyCam.end();
    
    // FFT //
	ofSetColor(255);
	ofPushMatrix();
	ofTranslate(16, 16);
    vid.draw(40,40);
	ofDrawBitmapString("Time Domain", 0, 0);
	
	soundMutex.lock();
	drawBuffer = middleBuffer;
	drawBins = middleBins;
	soundMutex.unlock();
	
	plot(drawBuffer, plotHeight / 2, 0);
	ofTranslate(0, plotHeight + 16);
	ofDrawBitmapString("Frequency Domain", 0, 0);
	plot(drawBins, -plotHeight, plotHeight / 2);
	ofTranslate(0, plotHeight + 16);
	spectrogram.update();
	spectrogram.draw(0, 0);
	ofDrawRectangle(0, 0, bufferSize, bufferSize / 2);
	ofDrawBitmapString("Spectrogram", 0, 0);
	ofPopMatrix();
	string msg = ofToString((int) ofGetFrameRate()) + " fps";
	ofDrawBitmapString(msg, appWidth - 80, appHeight - 20);
    
    // OSC //
    // draw recent unrecognized messages
    for(int i = 0; i < NUM_MSG_STRINGS; i++){
        ofDrawBitmapStringHighlight(msgStrings[i], 10, 40 + 15 * i);
    }
    
    string buf = "listening for osc messages on port " + ofToString(PORT);
    ofDrawBitmapStringHighlight(buf, 10, 20);
    
    // draw mouse state
    ofPoint mouseIn(mouseXf*ofGetWidth(), mouseYf*ofGetHeight());
    if(mouseButtonInt == 0){
        ofSetColor(255);
    } else{
        ofSetColor(ofColor::salmon);
    }
    ofDrawCircle(mouseIn, 20);
    ofDrawBitmapStringHighlight(mouseButtonState, mouseIn);
}


void ofApp::plot(vector<float>& buffer, float scale, float offset) {
	ofNoFill();
	int n = buffer.size();
	ofDrawRectangle(0, 0, n, plotHeight);
	glPushMatrix();
	glTranslatef(0, plotHeight / 2 + offset, 0);
	ofBeginShape();
	for (int i = 0; i < n; i++) {
		ofVertex(i, buffer[i] * scale);
	}
	ofEndShape();
	glPopMatrix();
}

void ofApp::audioReceived(float* input, int bufferSize, int nChannels) {
	if (mode == MIC) {
		// store input in audioInput buffer
		memcpy(&audioBuffer[0], input, sizeof(float) * bufferSize);
		
		float maxValue = 0;
		for(int i = 0; i < bufferSize; i++) {
			if(abs(audioBuffer[i]) > maxValue) {
				maxValue = abs(audioBuffer[i]);
			}
		}
		for(int i = 0; i < bufferSize; i++) {
			audioBuffer[i] /= maxValue;
		}
		
	} else if (mode == NOISE) {
		for (int i = 0; i < bufferSize; i++)
			audioBuffer[i] = ofRandom(-1, 1);
	} else if (mode == SINE) {
		for (int i = 0; i < bufferSize; i++)
			audioBuffer[i] = sinf(PI * i * mouseX / appWidth);
	}
	
	fft->setSignal(&audioBuffer[0]);

	float* curFft = fft->getAmplitude();
	memcpy(&audioBins[0], curFft, sizeof(float) * fft->getBinSize());

	float maxValue = 0;
	for(int i = 0; i < fft->getBinSize(); i++) {
		if(abs(audioBins[i]) > maxValue) {
			maxValue = abs(audioBins[i]);
		}
	}
	for(int i = 0; i < fft->getBinSize(); i++) {
		audioBins[i] /= maxValue;
	}
	
	int spectrogramWidth = (int) spectrogram.getWidth();
	int n = (int) spectrogram.getHeight();

	for(int i = 0; i < n; i++) {
		int j = (n - i - 1) * spectrogramWidth + spectrogramOffset;
		int logi = ofMap(powFreq(i), powFreq(0), powFreq(n), 0, n);
        spectrogram.setColor(j, (unsigned char) (255. * audioBins[logi]));
	}
	spectrogramOffset = (spectrogramOffset + 1) % spectrogramWidth;
	
	soundMutex.lock();
	middleBuffer = audioBuffer;
	middleBins = audioBins;
	soundMutex.unlock();
}

void ofApp::keyPressed(int key) {
	switch (key) {
	case 'm':
		mode = MIC;
		break;
	case 'n':
		mode = NOISE;
		break;
	case 's':
		mode = SINE;
		break;
    case 'f':
        frameByframe=!frameByframe;
        vid.setPaused(frameByframe);
        break;
    case OF_KEY_LEFT:
        vid.previousFrame();
        break;
    case OF_KEY_RIGHT:
        vid.nextFrame();
        break;
    case '0':
        vid.firstFrame();
        break;
	}
}


//--------------------------------------------------------------


//--------------------------------------------------------------
void ofApp::keyReleased(int key){
    
}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){
    if(!frameByframe){
        int width = ofGetWidth();
        float pct = (float)x / (float)width;
        float speed = (2 * pct - 1) * 5.0f;
        vid.setSpeed(speed);
    }
}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){
    if(!frameByframe){
        int width = ofGetWidth();
        float pct = (float)x / (float)width;
        vid.setPosition(pct);
    }
}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){
    if(!frameByframe){
        vid.setPaused(true);
    }
}


//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){
    if(!frameByframe){
        vid.setPaused(false);
    }
}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){
    
}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){
    
}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){
    
}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){
    
}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){
    
}
