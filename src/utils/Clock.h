/*
 * Clock.h
 *
 *  Created on: Jun 17, 2014
 *      Author: dsalinas
 */

#ifndef GUDHI_CLOCK_H_
#define GUDHI_CLOCK_H_

#include <sys/time.h>
#include <iostream>


class Clock{

public:
	Clock(){
		end_called = false;
		begin();
		msg = "time";
	}

	Clock(const std::string& msg_){
		end_called = false;
		begin();
		msg = msg_;
	}

	void begin() const{
		end_called = false;
		gettimeofday(&startTime, NULL);
	}

	void end() const{
		end_called = true;
		gettimeofday(&endTime, NULL);
	}

	void print() const{
		std::cout << *this << std::endl;
	}

	friend std::ostream& operator<< (std::ostream& stream,const Clock& clock){
		//		if(!clock.end_called) clock.end();
		if(!clock.end_called) clock.end();

		long totalTime =  (clock.endTime.tv_sec - clock.startTime.tv_sec) * 1000000L;
		totalTime += (clock.endTime.tv_usec - clock.startTime.tv_usec);
		stream << clock.msg <<":"<<clock.num_seconds() <<"s";

		return stream;
	}

	double num_seconds() const{
		if(!end_called) end();
		long totalTime =  (endTime.tv_sec - startTime.tv_sec) * 1000000L;
		totalTime += (endTime.tv_usec - startTime.tv_usec);
		return (totalTime / 1000L)/(double) 1000;
	}

private:
	mutable struct timeval startTime, endTime;
	mutable bool end_called;
	std::string msg;
};


#endif /* GUDHI_CLOCK_H_ */
