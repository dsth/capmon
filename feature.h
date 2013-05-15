#ifndef FEATURE_H
#define FEATURE_H

class feature_min { 

    unsigned int _gstart;
    unsigned int _gend;

  public:
    feature_min( unsigned int a,  unsigned int b) : _gstart(a), _gend(b) {} 
    virtual ~feature_min() {} // not sure if using derived classes via base pointer but better to avoid risk.
    unsigned int gstart() const { return _gstart; } 
    unsigned int gend() const { return _gend; } 

};

#endif
