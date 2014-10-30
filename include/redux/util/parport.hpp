#ifndef REDUX_UTIL_PARPORT_HPP
#define REDUX_UTIL_PARPORT_HPP

#include <cstdint>
#include <mutex>

// This maps parallel port pin numbers to register bit numbers.
// - shifts of 0..7   represent bits 0..7 of BASE + 0 (data register)
// - shifts of 8..15  represent bits 0..7 of BASE + 1 (status register)
// - shifts of 16..23 represent bits 0..7 of BASE + 2 (control register)

#define PAR_PIN01 (1 << 16)                     // ~STROBE          control_0       In/Out
#define PAR_PIN02 (1 << 0)                      // DATA             data_0-7        In/Out
#define PAR_PIN03 (1 << 1)
#define PAR_PIN04 (1 << 2)
#define PAR_PIN05 (1 << 3)
#define PAR_PIN06 (1 << 4)
#define PAR_PIN07 (1 << 5)
#define PAR_PIN08 (1 << 6)
#define PAR_PIN09 (1 << 7)
#define PAR_PIN10 (1 << 14)                     // ACK              status_6          In
#define PAR_PIN11 (1 << 15)                     // ~Busy            status_7          In
#define PAR_PIN12 (1 << 13)                     // Paper-out        status_5          In
#define PAR_PIN13 (1 << 12)                     // Select           status_4          In
#define PAR_PIN14 (1 << 17)                     // ~Line-feed       control_1       In/Out
#define PAR_PIN15 (1 << 11)                     // Error            status_3          In
#define PAR_PIN16 (1 << 18)                     // Reset            control_2       In/Out
#define PAR_PIN17 (1 << 19)                     // ~Select-Printer  control_3       In/Out
#undef  PAR_PIN18
#undef  PAR_PIN19 // pins 18..25 are commoned to signal ground
#undef  PAR_PIN20 // (not controllable)
#undef  PAR_PIN21
#undef  PAR_PIN22
#undef  PAR_PIN23
#undef  PAR_PIN24
#undef  PAR_PIN25

#define PAR_IRQ_MODE   (1 << 20) // control_4 ("pseudo-pin") -> pin 10 causes interrupts
#define PAR_INPUT_MODE (1 << 21) // control_5 ("pseudo-pin") -> data pins are input

#define PAR_MASK_PIN      0x0ff8ff      // Register-bits corresponding to physical pins
#define PAR_MASK_DATA     0x0000ff      // Register at baseaddress + 0 (data register)
#define PAR_MASK_STATUS   0x00ff00      // Register at baseaddress + 1 (status register) 
#define PAR_MASK_CONTROL  0xff0000      // Register at baseaddress + 2 (control register)
#define PAR_MASK_INVERT   0x0b8000      // Pins with inverted logic, i.e. register-low <=> high pin-voltage

 
namespace redux {

    namespace util {


        /*!  @ingroup util
         *  @{
         */


        /*!  @file      parport.hpp
         *   @brief     Convenience class for accessing a parallel port.
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2014
         */


        class ParPort {
            
        public:
            enum status_t { PAR_ERR=-1, PAR_NONE, PAR_OK };
            union regtype {
                uint32_t i;
                uint8_t b[4];
            };
            
            static uint32_t PIN[26];
            
            ParPort(uint64_t address=0);
            ~ParPort();
            
            void init(void);
            void release(void);
            void setMasks(void);
            void setBaseAddress(uint64_t);
            
            uint8_t read_register(uint8_t regnum);
            uint32_t read_registers(void);
            void write_register(uint8_t regnum, uint8_t value);
            void write_registers(uint32_t value);
            
            void clear_data(uint8_t value);
            void set_data(uint8_t value);
            void clear_control(uint8_t value);
            void set_control(uint8_t value);
            
            void print(void);
            
        protected:
            uint64_t baseaddress;
            regtype registers;   // b[0]=data, b[1]=status, b[2]=control
            status_t status;

            constexpr static regtype inv_mask {PAR_MASK_INVERT}; 
            regtype input_mask;
            regtype output_mask;
            
        private:
            std::mutex mutex;
           
        };
        
        /*! @} */

    }

}

#endif // REDUX_UTIL_PARPORT_HPP
