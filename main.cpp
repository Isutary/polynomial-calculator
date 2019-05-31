#include<iostream>
#include<cmath>
#include<algorithm>
#include "Polinom.h"
using namespace std;

// racuna presjek polinoma binarnom pretragom
double PresjekPolinoma(const Polinom &p1, const Polinom &p2, double a, double b, double e=1e-7) {
    if ((p1(a)-p2(a))*(p1(b)-p2(b)) > 0) throw invalid_argument("Polinomi se ne sijeku na intervalu [a,b] ili se sijeku vise puta!");
    double mid = (a+b)/2;
    if (p1(a) < p2(a)) {
        while(abs(p1(mid)-p2(mid)) > e) {
            if (p1(mid) < p2(mid)) a = mid;
            else b = mid;
            mid = (a+b)/2;
        }
    }
    else {
        while(abs(p1(mid)-p2(mid)) > e) {
            if (p1(mid) > p2(mid)) a = mid;
            else b = mid;
            mid = (a+b)/2;
        }
    }
    return mid;
}

int main() {

    vector<double> v1 {1,1,3,-5};
    vector<double> v2 {1,1,3,5};
    try {
        Polinom p1(v1);
        Polinom p2(v2);
        Polinom p3(3,2);
        Polinom p4;
        Polinom niz[]{ p1,p2,p3 };
        sort(niz, niz+3);
        cout << "Niz sortiranih polinoma je: " << endl;
        for (auto i : niz) cout << i << endl;
        cout << "Koeficijenti polinoma p1 su: ";
        for (unsigned int i=0; i<p1.size(); i++) cout << p1[i] << " ";
        cout << endl;
        cout << "Vrijednost polinoma p2 u tacki 2 je: " << p2(2) << endl;
        cout << "Kvadrat polinoma p3 je: " << (p3^2) << endl;
        cout << "Nula polinoma p2 je: " << p2.NulaPolinoma(-2, 100) << endl;
        cout << "Zbir polinoma p1 i p3 je: " << p1+p3 << endl;
        cout << "Razlika polinoma p3 i p2 je: " << p3-p2 << endl;
        cout << "Proizvod polinoma p1 i p2 je: " << p1*p2 << endl;
        cout << "Manji od polinoma p1 i p2 je: " << ((p1<p2)?"p1":(p2<p1)?"p2":"p1 == p2") << endl;
        cout << "Veci od polinoma p1 i p2 je: " << ((p1>p2)?"p1":(p2>p1)?"p2":"p1 == p2") << endl;
        cout << "Polinomi p2 i p3 su: " << ((p2 == p3)?"jednaki":"razliciti") << endl;
        cout << "Unesite polinom p4: ";
        cin >> p4;
        cout << "Unesite tacke a i b: ";
        double a, b;
        cin >> a >> b;
        cout << "Vas polinom je: " << p4 << endl;
        cout << "Vas polinom sjece polinom p1 u tacki: " << PresjekPolinoma(p1, p4, a, b) << endl;
    }
    catch (const invalid_argument &invArg) {
        cout << invArg.what();
    }
    catch (const logic_error &logErr) {
        cout << logErr.what();
    }

    return 0;
}
