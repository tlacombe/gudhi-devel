#ifndef __TEST_H
#define __TEST_H

#include <list>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>

using namespace std;
/**
 * Class to perform test
 */

class Test
{
private :
	string name;
	bool (*test)();

	string separation() const{
		return "+++++++++++++++++++++++++++++++++++++++++++++++++\n";
	}

	string print_between_plus(string& s) const{
		stringstream res;
		res << "+++++++++++++++++"<<s<<"+++++++++++++++++\n";
		return res.str();
	}


public:
	Test(string name_,bool (*test_)()){
		name=name_;
		test =test_;
	}

	bool run(){
		cout << print_between_plus(name);
		return test();
	}
	string getName(){
		return name;
	}
};


class Tests
{
private:
	list<Test> tests;

public:
	void add(string name_,bool (*test_)()){
		Test test(name_,test_);
		tests.push_back(test);
	}
	bool run(){
		bool tests_succesful(true);
		vector<bool> res;
		for (Test test : tests){
			res.push_back(test.run());
		}
		cout << "\n\n results of tests : "<<endl;
		int i=0;
		for (Test t : tests){
			cout << "Test "<<i<< " \""<<t.getName()<<"\"  -->  ";
			if (res[i++]) cout << "OK"<<endl;
			else {
				cout << "Fail"<<endl;
				tests_succesful = false;
				break;
			}
		}
		return tests_succesful;

	}
};

#endif
