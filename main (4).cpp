
#include "Matrix.h"

template <typename T>
Matrix<T>::Matrix()
{
    m = 0;
    n = 0;
}

template <typename T>
void Matrix<T>::Create(int i, int j)
{   try
    {
        if ((i <= 0) || (j <= 0)) throw MyExeption();
        else
        {
            m = i;
            n = j;
            a.reserve(m*n);
        }
    }
    catch(MyExeption)
    {
        cerr<<"Rows or columns must be > 0, matrix wasn't created."<<endl;
        return;
    }
}

template <typename T>
Matrix<T>::Matrix(int m, int n)
{
    Create(m, n);
}

template <typename T>
Matrix<T>::Matrix(int m, int n, vector<T> vec)
{
    Create_with_data(m , n, vec);
}





template <typename T>
void Matrix<T>::Create_with_data(int i, int j, vector<T> vec)
{
    try {
        if ((i <= 0) || (j <= 0)) throw (MyExeption());
        else {
            m = i;
            n = j;
            a.resize(vec.size());
            copy(vec.begin(), vec.end(), a.begin());
        }
    }
    catch(MyExeption)
    {
        cerr<<"Rows or columns must be > 0, matrix wasn't created."<<endl;
        return;
    }
}

template <typename T>
[[nodiscard]] pair<int, int> Matrix<T>::M_Size() const
{
    return make_pair(m, n);
}

template <typename T>
Matrix<T> Matrix<T>::AdamarMultiple(const Matrix<T> &obj)
{
    try {
        if ((this->M_Size().first != obj.M_Size().first) || (this->M_Size().second != obj.M_Size().second)) throw MyExeption();
        else {
            Matrix temp(this->M_Size().first, this->M_Size().second);
            for (int i = 0; i < this->M_Size().first; i++) {
                for (int j = 0; j < this->M_Size().second; j++) {
                    temp.a.push_back(this->a[n * i + j] * obj.a[n * i + j]);
                }
            }
            return temp;
        }
    }
    catch (MyExeption)
    {
        cerr << "Dimentions of matrix don't match, matrix wasn't changed." << '\n';
        return *this;
    }
}


template <typename T>
double Matrix<T>::Det()
{
    try {
        if (this->Rank() != this->M_Size().first || this->M_Size().first != this->M_Size().second)
            throw MyExeption();
        else
        {
            m = this->M_Size().first;
            n = this->M_Size().second;
            vector<T> copy = this->a;
            for (int k = 0; k < m - 1; k++)
            {
                if (this->a[k*n + k] == 0)
                {
                    int trigger = 0;
                    for (int r = k; r < m; r++)
                        if (this->a[r*n + k] != 0)
                        {
                            trigger = 1;
                            for (int j = k; j < m; j++)
                            {
                                this->a[k*n + j] += this->a[r*n + j];
                            }
                            break;
                        }
                    if (trigger == 0)
                        return 0;
                }
                double tmp_up = this->a[k*n + k];
                for (int i = k+1; i < m; i++)
                {
                    double tmp_low = this->a[i*n + k];
                    for (int j = k; j < m; j++)
                    {
                        this->a[i*n + j] += (this->a[k*n + j] * (-tmp_low)) / tmp_up;
                    }
                }
            }
            double det = 1;
            for (int d = 0; d < m; d++)
                det *= this->a[d*n + d];
            this->a = copy;
            return det;
        }
    }
    catch(MyExeption)
    {
        cerr << " This is not a square matrix or rank of matrix less than dimention, rank wasn't calculated" <<endl;
        return 0;
    }
}


template <typename T>
Matrix<T> Matrix<T>:: Inverse()
{
    try {
        if (this->Rank() != this->M_Size().first || this->M_Size().first != this->M_Size().second)
            throw MyExeption();
        else {
            int N = this->n;
            Matrix<T> minor;
            vector<T> vec1((N-1)*(N-1));
            minor.Create_with_data(N-1, N-1, vec1);
            Matrix<T> souz;
            vector<T> vec(N*N);
            souz.Create_with_data(N,N,vec);
            this->Transpose();
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++) {
                    for (int a = 0; a < N - 1; a++)
                        for (int b = 0; b < N - 1; b++) {
                            int p1 = a, p2 = b;
                            if (a >= i) p1 += 1;
                            if (b >= j) p2 += 1;
                            minor.a[a * (N - 1) + b] = this->a[p1 * N + p2];
                        }
                    souz.a[i * N + j] = minor.Det() * pow((-1), i + j);
                }
            souz * (1 / this->Det());
            this->Transpose();
            return souz;
        }
    }
    catch (MyExeption)
    {
        cerr << "This matrix doesn't have inverse matrix, matrix wasn't changed." << '\n';
        return *this;
    }
}

template <typename T>
int Matrix<T>::Rank()
{
    int i, j, k, l;
    double r;
    m = this->m;
    n = this->n;
    double eps = 0.000001;
    vector<T> a = this->a;
    i = 0; j = 0;
    while (i < m && j < n) {
        r = 0.0;
        for (k = i; k < m; ++k) {
            if (fabs(a[k*n + j]) > r) {
                l = k;
                r = fabs(a[k*n + j]);
            }
        }
        if (r <= eps) {
            for (k = i; k < m; ++k) {
                a[k*n + j] = 0.0;
            }
            ++j;
            continue;
        }

        if (l != i) {
            for (k = j; k < n; ++k) {
                r = a[i*n + k];
                a[i*n + k] = a[l*n + k];
                a[l*n + k] = (-r);
            }
        }
        for (k = i+1; k < m; ++k) {
            r = (-a[k*n + j] / a[i*n + j]);


            a[k*n + j] = 0.0;
            for (l = j+1; l < n; ++l) {
                a[k*n + l] += r * a[i*n + l];
            }
        }

        ++i; ++j;
    }
    return i;

}

template <typename T>
T Matrix<T>::Trace() const
{
    try
    {
        if (this->M_Size().first != this->M_Size().second) throw MyExeption();
        else {
            int r = this->M_Size().first;
            int c = this->M_Size().second;
            T trace = 0;
            for (int i = 0; i < r; i++) {
                for (int j = 0; j < c; j++) {
                    if (i == j ) trace += a[i * c + j];
                }
            }
            return trace;
        }
    }
    catch(MyExeption)
    {
        cerr <<" It is not a square matrix, trace wasn't calculated." << '\n';
        return NULL;
    }
}

template <typename T>
double Matrix<T>::FrobeniusNorm() const
{
    double norm = 0.0;
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            norm += pow(a[i * n+ j], 2);
        }
    }
    return pow(norm, 0.5);
}



template <typename T>
T Matrix<T>::ScalarMul(const Matrix<T> &obj) const
{
    try
    {
        if ((this->M_Size().first==1 || this->M_Size().second ==1) &&
            (obj.M_Size().first == 1 || obj.M_Size().second ==1)
            &&  (this->M_Size().first * this->M_Size().second == obj.M_Size().first * obj.M_Size().second)) {
            T scal = 0;
            int len = obj.M_Size().first * obj.M_Size().second;
            for (int i = 0; i < len; i++)
                scal += this->a[i] * obj.a[i];
            return scal;
        }
        else throw MyExeption();
    }
    catch (MyExeption)
    {
        cerr << "This is not a vector m*1 or 1*n, or size of vectors are not the same, scalar mul wasn't calculated."<<endl;
        return 0;
    }

}

template <typename T>
double Matrix<T>::EuclidNorm() const
{
    try
    {
        if ((this->M_Size().first == 1 || this->M_Size().second == 1)) {

            double norm = 0;
            int len = this->M_Size().first * this->M_Size().second;
            for (int i = 0; i < len; i++)
                norm += a[i] * a[i];
            norm = pow(norm, 0.5);
            return norm;
        }
        else throw MyExeption();

    }
    catch (MyExeption)
    {
        cerr<< " This is not a vector m*1 or 1*n, euclid norm wasn't calculated." << '\n';
        return 0;
    }
}
template <typename T>
T Matrix<T>::MaxNorm() const
{
    try {
        if (this->M_Size().first == 1 || this->M_Size().second == 1){

            T max = 0;
            for (int i = 0; i < this->M_Size().first * this->M_Size().second; i++)
                if (abs(this->a[i]) > max)
                    max = abs(a[i]);
            return max;
        }
        else throw MyExeption();

    }
    catch (MyExeption)
    {
        cerr<< "This is not a vector m*1 or 1*n, max norm wasn't calculated."<< '\n';
        return 0;
    }
}

template <typename T>
Matrix<T> Matrix<T>::VectorMul(const Matrix<T> &obj)
{
    try
    {
        if ((this->M_Size().first==1 || this->M_Size().second ==1) &&
            (obj.M_Size().first == 1 || obj.M_Size().second ==1) &&
            (this->M_Size().first * this->M_Size().second == obj.M_Size().first * obj.M_Size().second) &&
            (this->M_Size().first==3 || this->M_Size().second ==3) &&
            (obj.M_Size().first == 3 || obj.M_Size().second ==3))
        {
            Matrix<T> c(1, 3);
            c.a.push_back((this->a[1] * obj.a[2]) - (obj.a[1] * this->a[2]));
            c.a.push_back(-(this->a[0] * obj.a[2]) + (obj.a[0] * this->a[2]));
            c.a.push_back((this->a[0] * obj.a[1]) - (obj.a[0] * this->a[1]));
            return c;
        }
        else throw MyExeption();
    }
    catch (MyExeption)
    {
        cerr << "This is not a vector m*1 or 1*n, or size of vectors are not the 3D, vector mul wasn't calculated."<< '\n';
        return *this;
    }

}
template <typename T>
double Matrix<T>::AngelofVec(const Matrix<T> &obj)
{
    try
    {
        if ((this->M_Size().first==1 || this->M_Size().second ==1) &&
            (obj.M_Size().first == 1 || obj.M_Size().second ==1)
            &&  (this->M_Size().first * this->M_Size().second == obj.M_Size().first * obj.M_Size().second)) {
            double cos;
            cos = (this->ScalarMul(obj)) / (this->EuclidNorm() * obj.EuclidNorm());
            return (acos(cos) * 180.0) / 3.14159265;
        }
        else throw MyExeption();
    }
    catch (MyExeption)
    {
        cerr<< "These are not vectors m*1 or 1*n, angel wasn't calculated"<<'\n';
        return 0;
    }
}

template <typename T>
void Matrix<T>::Transpose()
{
    vector<T> vec = this->a;
    int t = this->m;
    this->m = this->n;
    this->n = t;
    for (int i = 0; i < this->m; i++)
        for (int j = 0; j < this->n; j++)
            vec[i* this->n + j] = this->a[j * this->n + i];
    this->a.clear();
    this->a = vec;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &obj) {

    try {
        if (this->M_Size().second != obj.M_Size().first) throw MyExeption();
        else {
            T el = 0;
            Matrix temp(this->M_Size().first, obj.M_Size().second);
            for (int i = 0; i < this->M_Size().first; i++) {
                for (int j = 0; j < obj.M_Size().second; j++) {
                    for (int k = 0; k < this->M_Size().second; k++) {
                        el += (this->a[i * this->M_Size().second + k] * obj.a[k * obj.M_Size().second + j]);
                    }
                    temp.a.push_back(el);
                    el = 0;

                }
            }
            return temp;
        }
    }
    catch (MyExeption) {
        cerr << "Dimentions of matrix don't match( m * n x n * l), matrix wasn't changed." << '\n';
        return *this;
    }


}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &obj) {
    try {
        if ((this->M_Size().first != obj.M_Size().first) || (this->M_Size().second != obj.M_Size().second)) throw MyExeption();
        else {
            Matrix temp(this->M_Size().first, this->M_Size().second);
            for (int i = 0; i < this->M_Size().first; i++) {
                for (int j = 0; j < this->M_Size().second; j++) {
                    temp.a.push_back(this->a[this->M_Size().second * i + j] + obj.a[this->M_Size().second * i + j]);
                }
            }
            return temp;
        }
    }
    catch (MyExeption) {
        cerr << "Dimentions of matrix don't match, matrix wasn't changed."<<endl;
        return *this;
    }
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &obj) {
    try {
        if ((this->M_Size().first != obj.M_Size().first) || (this->M_Size().second != obj.M_Size().second)) throw MyExeption();
        else {
            Matrix temp(this->M_Size().first, this->M_Size().second);
            for (int i = 0; i < this->M_Size().first; i++) {
                for (int j = 0; j < this->M_Size().second; j++) {
                    temp.a.push_back(this->a[this->M_Size().second * i + j] - obj.a[this->M_Size().second * i + j]);
                }
            }
            return temp;
        }
    }
    catch (MyExeption) {
        cerr << "Dimentions of matrix don't match, matrix wasn't changed."<<endl;
        return *this;
    }
}

template <typename T>
Matrix<T> operator*(Matrix<T> &obj, T r)
{
    for (int i = 0; i < obj.M_Size().first * obj.M_Size().second; i++) {
        obj.a[i] *=  r;
    }
    return obj;
}

template <typename T>
Matrix<T> Matrix<T>::operator/(T r)
{
    for (int i = 0; i < this->M_Size().first * this->M_Size().second; i++) {
        this->a[i] /=  r;
    }
    return *this;
}

template <typename T>
Matrix<T> operator*(T r, Matrix<T> &obj)
{
    return obj * r;
}


template <typename T>
Matrix<T>&  Matrix<T>::operator=(const Matrix<T> &obj)
{
    if (&obj == this)  return *this;
    this->m = obj.M_Size().first;
    this->n = obj.M_Size().second;
    this->a.resize(obj.M_Size().first * obj.M_Size().second );
    for (int i = 0; i < obj.M_Size().first; i++)
        for (int j = 0; j < obj.M_Size().second; j++)
            this->a[i*obj.M_Size().second + j] = obj.a[i*obj.M_Size().second + j];

    return *this;
}



template <typename T>
ostream& operator<<(ostream& out, const Matrix<T>& matrix)
{
    for (int row = 0; row < matrix.M_Size().first; row++)
    {
        cout << "|";
        for (int column = 0; column < matrix.M_Size().second; column++)
        {
            if (column > 0)
            {
                out << " ";
            }
            out << matrix.a[row*matrix.M_Size().second + column];
        }
        cout << "|\n";
    }
    return out;
}

template <typename T>
istream&  operator>>(istream& in, Matrix<T>& matrix)
{
    T b;
    for (int i = 0; i < matrix.M_Size().first; i++)
        for (int j = 0; j < matrix.M_Size().second; j++) {
            in >> b;
            matrix.a.push_back(b);
        }
    return in;
}

template <typename T>
void Matrix<T>::write_to_bin_file(ostream& stream) const
{
    int n = this->n;
    int m = this->m;
    stream.write(reinterpret_cast<char*>(&n), sizeof(int));
    stream.write(reinterpret_cast<char*>(&m), sizeof(int));
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            T value = a[i*m + j];
            stream.write(reinterpret_cast<char*>(&value), sizeof(T));
        }
    }
}
template <typename T>
void Matrix<T>::read_from_bin_file(istream& stream)
{
    try {
        int r, c;
        stream.read(reinterpret_cast<char *>(&r), sizeof(int));
        stream.read(reinterpret_cast<char *>(&c), sizeof(int));
        if (r <= 0 || c <= 0) throw MyExeption();
        vector<T> vec(r * c);
        this->Create_with_data(r, c, vec);
        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < c; ++j) {
                stream.read(reinterpret_cast<char *>(&this->a[i * c + j]), sizeof(T));
                if (!stream.good()) throw MyExeption();
            }
        }
    }
    catch (MyExeption)
    {
        cerr << "Error of file reading." << endl;
        return;
    }
}
template <typename T>
void Matrix<T>::operator>> (ofstream &file) const
{
    if (this->m <= 0 || this->n <= 0) throw MyExeption();
    file << this->m << '\n' << this->n << '\n';
    for (int i = 0; i < this->m; ++i)
    {
        for (int j = 0; j < this->n; ++j)
        {
            file << this->a[i*this->n + j] << ' ';
        }
        file << '\n';
    }
}

template <typename T>
void Matrix<T>::operator<< (ifstream &file)
{
    int r, c;
    file >> r >> c;
    if (r <= 0 || c <= 0) throw MyExeption();
    vector<T> vec(r * c);
    this->Create_with_data(r, c, vec);
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {

            file >> this->a[i * c +j];
        }
    }

}

template <typename T>
void PCA<T>::Centering()
{
    int numr = this->mat.M_Size().first;
    int numc = this->mat.M_Size().second;
    double mean;
    for (int i =0 ; i < numc; ++i)
    {
        mean = 0.0;
        for(int j =0; j < numr; ++j)
        {
            mean += this->mat.a[j*numc + i];
        }
        mean = mean/numr;
        for(int k =0; k < numr; ++k)
        {
            this->mat.a[k*numc + i] -= mean;
        }
    }
}

template <typename T>
void PCA<T>::Normalize()
{
    int numr = this->mat.M_Size().first;
    int numc = this->mat.M_Size().second;
    double mean;
    for (int i =0 ; i < numc; ++i)
    {
        mean = 0.0;
        for(int j =0; j < numr; ++j)
        {
            mean += this->mat.a[j*numc + i];
        }
        mean = mean/numr;
        double diviation = 0.0;
        for(int k =0; k < numr; ++k)
        {
            diviation += pow((this->mat.a[k*numc + i] - mean),2);
        }
        diviation /= numr;
        diviation = pow(diviation, 0.5);
        for(int k =0; k < numr; ++k)
        {
            this->mat.a[k*numc + i] = (this->mat.a[k*numc + i] - mean)/diviation;
        }
    }
}
template <typename T>
void Matrix<T>::Transpose_vec()
{
    int temp = this->m;
    this->m = this->n;
    this->n = temp;

}

template <typename T>
vector<Matrix<T>> PCA<T>::NIPALS(int num_PC)
{
    vector<Matrix<T>> result;
    this->Centering();
    this->Normalize();
    int numr = this->mat.M_Size().first;//
    int numc = this->mat.M_Size().second;//
    int PC = num_PC;

    Matrix<T> copy = this->mat;
    Matrix<T> scores;
    Matrix<T> loadings;

    Matrix<T> told(1, numr); //1 * r
    Matrix<T> p;
    Matrix<T> tnew;

    vector<T> v(numr * PC);
    vector<T> q(numc * PC);
    vector<T> yo(numr);

    tnew.Create_with_data(1,numr, yo);
    scores.Create_with_data(numr, PC, v); // r * c
    loadings.Create_with_data(numc, PC, q);

    double eps = 0.00000001;
    int ij = 0;
    int flag = 0;

    for (int i = 0; i < PC; ++i) //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    {
        while (1)
        {
            if (flag == 0)
            {
                for (int j = 0; j < numr; ++j)
                {
                    told.a.push_back(copy.a[j*numc + i]);
                }
            }
            tnew.a = told.a;
            p = tnew * copy / tnew.ScalarMul(tnew); //1 * r * r * c = 1 * c
            p / p.EuclidNorm();// 1 * c
            p.Transpose_vec();// c * 1
            tnew = copy * p / p.ScalarMul(p);//r * c *  c * 1 = r * 1
            tnew.Transpose_vec();// 1 * r
            Matrix<T> d = told - tnew;
            if (d.EuclidNorm() <= eps) {
                tnew.Transpose_vec();//r * 1
                p.Transpose_vec();//1*c
                copy = copy - tnew * p;
                told.a.clear();
                flag = 0;
                break;
            }
            else
            {
                flag = 1;
                told.a = tnew.a;

            }
        }
        for (int j = 0; j < numr; j++) {
            scores.a[j*PC + ij] = tnew.a[j];
        }
        for (int j = 0; j < numc; j++) {
            loadings.a[j*PC + ij] = p.a[j];
        }
        ij += 1;
        tnew.Transpose_vec();
    }
    scores * (-1.0);
    loadings * (-1.0);
    result.push_back(scores);
    result.push_back(loadings);
    result.push_back(copy);
    return result ;

}
template <typename T>
double PCA<T>::TRV()
{
    int numr = this->m;
    int numc = this->n;
    double sum = 0;
    for(int i = 0; i < numr; i++)
    {
        for(int j = 0;j < numc; j++)
        {
            sum += pow(this->a[i*numc+j], 2);
        }

    }
    sum /= numr;
    sum /= numc;
    return sum;
}
template <typename T>
double PCA<T>::ERV()
{
    return (1 - this->TRV());
}
template <typename T>
double PCA<T>::Leverage(int num)
{
    Matrix<T> copy = *this;
    int numr = this->m;
    int numc = this->n;
    vector<T> t(numc);
    Matrix<T> tm(numc,1,t);
    for(int i =0;i < numc; i ++)
    {
        tm.a[i] = (this->a[num * numc + i]);
    }
    Matrix<T> tn = tm;
    tm.Transpose_vec();
    copy.Transpose();
    Matrix<T> pip = (copy * (*this));
    tm = tm * pip.Inverse() * tn;
    return tm.a[0];
}


template <typename T>
UnitMatrix<T>::UnitMatrix(int dim)
{
    Matrix<T>::Create(dim, dim);
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            if (i == j) this->a.push_back(1);
            else this->a.push_back(0);
}

template <typename T>
DiagonalMatrix<T>::DiagonalMatrix(int dim, vector<T> vec)
{
    Matrix<T>::Create(dim, dim);
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            if (i == j) this->a.push_back(vec[i]);
            else this->a.push_back(0);
}

int main()
{/*Matrix<double> x(2,2,{6.9, 9.6 ,1.1 ,4.4});
    Matrix<double> y(2,2);
    cin >> y;
    fstream file;
    file.open("matBinData.bin", ios::out | ios::binary);
    x.write_to_bin_file(file);
    y.write_to_bin_file(file);
    file.close();
    Matrix<double> r;
    Matrix<double> t;
    file.open("matBinData.bin", ios::in | ios::binary);
    r.read_from_bin_file(file);
    std::cout << r;
    t.read_from_bin_file(file);
    std::cout << t;*/
    /*Matrix<double> m(7,4);
    cin >> m;
    ofstream file("matrix1.txt");
    m >> file;
    file.close();
    ifstream file1("matrix1.txt");
    Matrix<double> t;
    t << file1;
    cout << t;
    file1.close();*/
    ifstream file1("data.txt");
    Matrix<double> t;
    t << file1;
    file1.close();
    PCA<double> t_(t);
    vector<Matrix<double>> result = t_.NIPALS(10);
    cout << result[0];
    //result[1].Transpose();
    //cout << result[0] * result[1] + result[2] ;
    //cout << result[0] <<'\n'<< result[1] << '\n' <<  result[2] ;
    //cout << result[0].Leverage(2);


}
