#ifndef POLINOM_CPP_INCLUDED
#define POLINOM_CPP_INCLUDED
#include "Polinom.h"
#include<iostream>
#include<regex>
#include<cmath>
using namespace std;

void Polinom::inputFun(string s1) {
    s1.erase(remove(s1.begin(), s1.end(), ' '), s1.end());
    skup.clear();
    for (unsigned int i=0; i<s1.size(); i++) {
        unsigned int step=0;
        double koef=0;
        string monom;
        int sign = min(s1.find('-', i+1), s1.find('+', i+1));
        if (sign == -1) {
            monom = s1.substr(i);
            i = s1.size()-1;
        }
        else {
            monom = s1.substr(i, sign-i);
            i = sign-1;
        }
        if (monom[0] != '+' && monom[0] != '-') monom = "+"+monom;
        if (regex_match(monom, regex("\\+(\\d+)x\\^(\\d+)")) || regex_match(monom, regex("-(\\d+)x\\^(\\d+)"))) {
            string::size_type n;
            koef = stod(monom, &n);
            monom = monom.substr(n+2);
            step = stoi(monom);
        }
        else if (regex_match(monom, regex("\\+(\\d+)\\*x\\^(\\d+)")) || regex_match(monom, regex("-(\\d+)\\*x\\^(\\d+)"))) {
            string::size_type n;
            koef = stod(monom, &n);
            monom = monom.substr(n+3);
            step = stoi(monom);
        }
        else if (regex_match(monom, regex("\\+x\\^(\\d+)"))) {
            koef = 1;
            monom = monom.substr(3);
            step = stoi(monom);
        }
        else if (regex_match(monom, regex("-x\\^(\\d+)"))) {
            koef = -1;
            monom = monom.substr(3);
            step = stoi(monom);
        }
        else if (regex_match(monom, regex("\\+x"))) {
            koef = 1;
            step = 1;
        }
        else if (regex_match(monom, regex("-x"))) {
            koef = -1;
            step = 1;
        }
        else if (regex_match(monom, regex("\\+(\\d+)x")) || regex_match(monom, regex("-(\\d+)x"))) {
            koef = stod(monom);
            step = 1;
        }
        else if (regex_match(monom, regex("\\+(\\d+)\\*x")) || regex_match(monom, regex("-(\\d+)\\*x"))) {
            koef = stod(monom);
            step = 1;
        }
        else if (regex_match(monom, regex("\\+(\\d+)")) || regex_match(monom, regex("-(\\d+)"))) {
            koef = stod(monom);
            step = 0;
        }
        else {
            skup.clear();
            throw logic_error("Pogresno unesen polinom!");
        }
        if (step < skup.size()) skup[step] += koef;
        else {
            skup.resize(step+1);
            skup[step] += koef;
        }
    }
    checkSkup();
    nonZero();
}

void Polinom::checkSkup() {
    int newSize = skup.size();
    while (skup[newSize-1] == 0 && newSize > 0) newSize--;
    skup.resize(newSize);
}

void Polinom::nonZero() const {
    monomCount = 0;
    for (unsigned int i=0; i<skup.size(); i++) if (skup[i] != 0) monomCount++;
}

double Polinom::operator()(double x) const {
    double sum = 0;
    for (unsigned int i=0; i<size(); i++) sum += skup[i]*pow(x,i);
    return sum;
}

Polinom Polinom::operator^(int n) const {
    if (n < 0) throw invalid_argument("Stepen ne moze biti manji od 0!");
    if (n == 0) return Polinom(0);
    if (n == 1) return (*this);
    int mid = n/2;
    Polinom tempP = (*this)^mid;
    if (n % 2 == 0) return tempP*tempP;
    else return tempP*tempP*(*this);
}

double Polinom::NulaPolinoma(double a, double b, double e) const {
    if (((*this)(a) < 0 && (*this)(b) < 0) || ((*this)(a) > 0 && (*this)(b) > 0)) throw invalid_argument("Polinom nema nula na intervalu [a,b] ili ima vise nula!");
    double mid = (a+b)/2;;
    if ((*this)(a) <= 0) {
        while (abs((*this)(mid)) > e) {
            if ((*this)(mid) < 0) a = mid;
            else b = mid;
            mid = (a+b)/2;
        }
    }
    else {
        while (abs((*this)(mid)) > e) {
            if ((*this)(mid) > 0) a = mid;
            else b = mid;
            mid = (a+b)/2;
        }
    }
    return mid;
}

Polinom operator+(const Polinom &p1, const Polinom &p2) {
    Polinom tempP;
    if (p1.size() < p2.size()) {
        tempP = p2;
        for (unsigned int i=0; i<p1.size(); i++) tempP[i] += p1[i];
    }
    else {
        tempP = p1;
        for (unsigned int i=0; i<p2.size(); i++) tempP[i] += p2[i];
    }
    tempP.checkSkup();
    tempP.nonZero();
    if (tempP.size() == 0) return Polinom();
    return tempP;
}

Polinom operator-(const Polinom &p1, const Polinom &p2) {
    Polinom tempP = p2;
    for (unsigned int i=0; i<tempP.size(); i++) tempP[i] = -tempP[i];
    return tempP + p1;
}

Polinom operator*(const Polinom &p1, const Polinom &p2) {
    if (p1.monomCount == 0 || p2.monomCount == 0) return Polinom();
    else if (p1.monomCount == 1) {
        Polinom tempP(p1.size()+p2.size()-2);
        for (unsigned int i = p1.size()-1; i < tempP.size(); i++) tempP[i] = p2[i - p1.size() + 1] * p1[p1.size() - 1];
        tempP.monomCount = p2.monomCount;
        return tempP;
    }
    else if (p2.monomCount == 1) {
        Polinom tempP(p1.size()+p2.size()-2);
        for (unsigned int i = p2.size()-1; i < tempP.size(); i++) tempP[i] = p1[i - p2.size() + 1] * p2[p2.size() - 1];
        tempP.monomCount = p1.monomCount;
        return tempP;
    }
    else {
        int midP1 = p1.size() / 2, midP2 = p2.size() / 2;
        int midP = min(midP1, midP2);
		Polinom A0(midP - 1), A1(p1.size() - midP - 1), B0(midP - 1), B1(p2.size() - midP - 1);
		for (unsigned int i = 0; i < A0.size(); i++) A0[i] = p1[i];
		A0.nonZero();
		for (unsigned int i = 0; i < A1.size(); i++) A1[i] = p1[i + midP];
		A1.nonZero();
		for (unsigned int i = 0; i < B0.size(); i++) B0[i] = p2[i];
		B0.nonZero();
		for (unsigned int i = 0; i < B1.size(); i++) B1[i] = p2[i + midP];
		B1.nonZero();
        Polinom Y = (A0+A1)*(B0+B1);
		Polinom U = A0*B0;
		Polinom Z = A1*B1;
		Polinom x(midP);
 		return U + (Y - U - Z)*x + Z*x*x;
    }
}

bool operator<(const Polinom &p1, const Polinom &p2) {
    if (p1[p1.size()-1] * p2[p2.size()-1] < 0) {
        if (p1[p1.size()-1] < 0) return true;
        else return false;
    }
    else if (p1.size() < p2.size()) return true;
    else if (p1.size() > p2.size()) return false;
    else {
        for (unsigned int i=0; i<p1.size(); i++) {
            if (p1[i] < p2[i]) return true;
            else if (p1[i] > p2[i]) return false;
        }
    }
    return false;
}

bool operator>(const Polinom &p1, const Polinom &p2) {
     if (p1[p1.size()-1] * p2[p2.size()-1] < 0) {
        if (p1[p1.size()-1] < 0) return false;
        else return true;
    }
    else if (p1.size() > p2.size()) return true;
    else if (p1.size() < p2.size()) return false;
    else {
        for (unsigned int i=0; i<p1.size(); i++) {
            if (p1[i] > p2[i]) return true;
            else if (p1[i] < p2[i]) return false;
        }
    }
    return false;
}

bool operator==(const Polinom &p1, const Polinom &p2) {
    return !(p1 < p2 || p1 > p2);
}

istream& operator>>(istream &flow, Polinom &p1) {
    string input;
    getline(cin, input);
    p1.inputFun(input);
    return flow;
}

ostream& operator<<(ostream &flow, const Polinom &p1) {
    if (p1.size() == 0) return flow;
    if (p1.size() == 1) { cout << p1[0]; return flow; }
    if (p1.size() == 2) {
        if (p1[1] == 1) cout << "x";
        else if (p1[1] == -1) cout << "-x";
        else cout << p1[1] << "x";
        if (p1[0] > 0) cout << "+" << p1[0];
        else if (p1[0] < 0) cout << p1[0];
        return flow;
    }
    if (p1[p1.size()-1] == 1) cout << "x^" << p1.size()-1;
    else if (p1[p1.size()-1] == -1) cout << "-x^" << p1.size()-1;
    else { cout << p1[p1.size()-1] << "x^" << p1.size()-1; }
    for(int i=p1.size()-2; i>=0; i--) {
        if (p1[i] != 0) {
            if (p1[i] == 1) {
                if (i == 1) cout << "+x";
                else if (i == 0) cout << "+1";
                else cout << "+x^" << i;
            }
            else if (p1[i] == -1) {
                if (i == 1) cout << "-x";
                else if (i == 0) cout << "-1";
                else cout << "-x^" << i;
            }
            else if (p1[i] > 0) {
                if (i == 1) cout << "+" << p1[i] << "x";
                else if (i == 0) cout << "+" << p1[i];
                else cout << "+" << p1[i] << "x^" << i;
            }
            else if (p1[i] < 0) {
                if (i == 1) cout << p1[i] << "x";
                else if (i == 0) cout << p1[i];
                else cout << p1[i] << "x^" << i;
            }
        }
    }
    return flow;
}

#endif // POLINOM_CPP_INCLUDED
