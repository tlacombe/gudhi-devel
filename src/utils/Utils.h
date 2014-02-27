/*
 * Utils.h
 *
 *  Created on: 13 juin 2013
 *      Author: salinasd
 */

#ifndef GUDHI_UTILS_H_
#define GUDHI_UTILS_H_


#define PRINT(a) std::cout << #a << ": " << (a) << " (DISP)"<<std::endl

#define DBG_VERBOSE
#ifdef DBG_VERBOSE
#define DBG(a) std::cout << "DBG: " << (a)<<std::endl
#define DBGMSG(a,b) std::cout << "DBG: " << a<<b<<std::endl
#define DBGVALUE(a) std::cout << "DBG: " <<  #a << ": " << a<<std::endl
#else
#define DBG(a)
#define DBGMSG(a,b)
#define DBGVALUE(a)
#endif




#endif /* UTILS_H_ */
