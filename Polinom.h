#ifndef POLINOM_H_INCLUDED
#define POLINOM_H_INCLUDED
#include<vector>
#include<stdexcept>
#include<string>
using namespace std;

class Polinom {
private:
    mutable int monomCount = 0; // broj nenultih koeficijenata u polinomu(vektoru)
    vector<double> skup; // vector koeficijenata polinoma
    void inputFun(string s1); // provjerava da li je polinom unesen u ispravnom obliku
    void checkSkup(); // smanji vektor, po potrebi, tako da je zadnji elemenat(koeficijent) nenulti
    void nonZero() const;  // broji nenulte elemenate(koeficijenate) u vektoru(polinomu)
public:
    Polinom() { skup.resize(1); skup[0]=0; } // konstruktor nula polinoma(nule)
    Polinom(int n, int k=1) : monomCount(1) { if (k == 0 || n < 0) throw invalid_argument("Neispravan stepen i/ili koeficijent!"); skup.resize(n+1); skup[n]=k; } // konstruktor za jednoclani polinom gdje je 'n' stepen, a 'k' koeficijent
    Polinom(vector<double> v) : skup(v) { checkSkup(); nonZero(); } // konstruktor preko vectora, uveden radi lakseg testiranja
    double& operator[](int i) { return skup[i]; } // uveden radi lakseg pristupa elementima polinoma
    const double& operator[](int i) const { return skup[i]; } // uveden radi lakseg pristupa elementima konstantnog polinoma
    unsigned int size() const { return skup.size(); } // uveden radi lakseg koristenja polinoma u petljama
    double operator()(double x) const; // racuna vrijednos polinoma u datoj tacki
    Polinom operator^(int n) const; // stepenuje polinom strategijom podijeli pa vladaj
    double NulaPolinoma(double a, double b, double e=1e-7) const; // racuna nulu polinoma sa slobodom 'e' binarnom pretragom
    friend Polinom operator+(const Polinom &p1, const Polinom &p2); // uveden radi lakse implementacije operatora '*'
    friend Polinom operator-(const Polinom &p1, const Polinom &p2); // uveden radi lakse implementacije operatora '*'
    friend Polinom operator*(const Polinom &p1, const Polinom &p2); // mnozi dva polinoma strategijom podijeli pa vlada
    friend bool operator<(const Polinom &p1, const Polinom &p2); // omogucava upotrebu funkcije 'sort'
    friend bool operator>(const Polinom &p1, const Polinom &p2); // uvedeen radi kompletnosti
    friend bool operator==(const Polinom &p1, const Polinom &p2); // uveden radi kompletnosti
    friend istream& operator>>(istream &flow, Polinom &p1); // omogucava unos polinoma operatorom '>>' pozivom funkcije 'inputFun'
    friend ostream& operator<<(ostream &flow, const Polinom &p1); // omogucava ispis polinoma operatorom '<<'
};

//#include "Polinom.cpp"
#endif // POLINOM_H_INCLUDED
