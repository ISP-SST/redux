#include "redux/util/parport.hpp"

#include <redux/util/stringutil.hpp>

#include <iostream>
#include <sys/io.h>

using namespace redux::util;
using namespace std;

#define INB(register_number) inb(baseaddress + register_number)
#define OUTB(value, register_number) outb(value, baseaddress + register_number)

// pins that are always inputs
#define PAR_ALWAYS_INPUT_PINS (PAR_PIN10 | PAR_PIN11 | PAR_PIN12 | PAR_PIN13 | PAR_PIN15)

// pins that can act as either inputs or outputs depending on PAR_INPUT_MODE
#define PAR_DATA_PINS (PAR_PIN02 | PAR_PIN03 | PAR_PIN04 | PAR_PIN05 | PAR_PIN06 | PAR_PIN07 | PAR_PIN08 | PAR_PIN09)

/* pins that can act as either input or output */
#define PAR_SWITCHABLE_PINS (PAR_PIN01 | PAR_PIN14 | PAR_PIN16 | PAR_PIN17)

uint32_t ParPort::PIN[] = {
        0,        // "pin 0" = undefined
        PAR_PIN01,
        PAR_PIN02,
        PAR_PIN03,
        PAR_PIN04,
        PAR_PIN05,
        PAR_PIN06,
        PAR_PIN07,
        PAR_PIN08,
        PAR_PIN09,
        PAR_PIN10,
        PAR_PIN11,
        PAR_PIN12,
        PAR_PIN13,
        PAR_PIN14,
        PAR_PIN15,
        PAR_PIN16,
        PAR_PIN17,
        0,0,0,0,0,0,0,0         // pin 18-25 - common ground, not accessible
    };

ParPort::ParPort ( uint64_t base ) : baseaddress ( base ), status(PAR_NONE) {
    if(baseaddress) init();
}


ParPort::~ParPort() {
    release();
}


void ParPort::init(void) {
    
    registers.i = 0;
    
    if ( !baseaddress ) {
        status = PAR_NONE;
        return;
    }
    
    
    if ( ioperm ( baseaddress,3,1 ) ) {
        status = PAR_ERR;
        return;
    }

    status = PAR_OK;


    read_registers();
    
    setMasks();
    
}


void ParPort::release(void) {
    
    if ( baseaddress && ioperm ( baseaddress,3,0 ) ) {
        status = PAR_ERR;
        return;
    }
    status = PAR_NONE;
}


void ParPort::setMasks(void) {
    output_mask.i = PAR_INPUT_MODE | PAR_IRQ_MODE;
    input_mask.i = PAR_ALWAYS_INPUT_PINS | PAR_INPUT_MODE | PAR_IRQ_MODE;
    if(registers.i&PAR_INPUT_MODE) input_mask.i |= PAR_DATA_PINS;
    else output_mask.i |= PAR_DATA_PINS;
}


void ParPort::setBaseAddress(uint64_t port) {
    release();
    baseaddress = port;
    init();
}


uint8_t ParPort::read_register ( uint8_t regnum ) {
    if( status != PAR_OK ) return 0;
    std::unique_lock<std::mutex> lock(mutex);
    registers.b[regnum] = INB ( regnum );
    //cout << "read_register(" << regnum << "," << bitString(registers.b[regnum]) << ")" << endl;
    return registers.b[regnum];
}


uint32_t ParPort::read_registers ( void ) {
    if( status != PAR_OK ) return 0;
    std::unique_lock<std::mutex> lock(mutex);
    for(uint8_t i = 0; i < 3; ++i) {
        registers.b[i] = INB ( i );
    }
    //cout << "read_registers(" << bitString(registers.i) << ")" << endl;
    return registers.i;
}


// store the written value in the lp_register variable, correct for
// inverted logic on some bits, and write a byte to an I/O port (base
// + regnum)
void ParPort::write_register ( uint8_t regnum, uint8_t value ) {
    if( status != PAR_OK ) return;
    std::unique_lock<std::mutex> lock(mutex);
#ifdef _DEBUG
    cout << "write_register(" << (int)regnum << "," << bitString(value) << ")" << endl;
#endif
    registers.b[regnum] = value;
    //if( regnum == 2 && value&(1<<5)) setMasks();
    OUTB ( value, regnum );
}


void ParPort::write_registers ( uint32_t value ) {
    if( status != PAR_OK ) return;
    std::unique_lock<std::mutex> lock(mutex);
#ifdef _DEBUG
    cout << "write_registers(" << bitString(value) << ")" << endl;
#endif
    registers.i = value;
    for(uint8_t i = 0; i < 3; ++i) {
        OUTB ( registers.b[i], i );
    }
}


void ParPort::clear_data(uint8_t value) {
    write_register(read_register(0) & ~value, 0);
}


void ParPort::set_data(uint8_t value) {
    write_register(read_register(0) | value, 0);
}


void ParPort::clear_control(uint8_t value) {
    
  // the control register requires:
  //    - switchable pins that are currently being used as inputs must be 1
  //    - all other pins may be either set or cleared

  uint8_t tmp = read_register(2);                               // read existing register
  tmp &= (~value);;                                             // set the requested bits to 1, leaving the others unchanged
  tmp |= (0x0f & ((input_mask.i & PAR_SWITCHABLE_PINS) >> 16));   // set all inputs to one (they may have been read as 0's!)

  write_register(tmp, 2);
  
}


void ParPort::set_control(uint8_t value) {
    
  // the control register requires:
  //    - switchable pins that are currently being used as inputs must be 1
  //    - all other pins may be either set or cleared

  uint8_t tmp = read_register(2);                               // read existing register
  tmp |= value;                                                 // set the requested bits to 1, leaving the others unchanged
  tmp |= (0x0f & ((input_mask.i & PAR_SWITCHABLE_PINS) >> 16));   // set all inputs to one (they may have been read as 0's!)

  write_register(tmp, 2);
  
}


void ParPort::print(void) {
    for(int i=0; i<3; ++i) cout << bitString(registers.b[i]) << endl;
}

