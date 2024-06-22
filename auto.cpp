#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <functional>
#include <cstdint>
#include <array>
#include <string>
#include <fstream>
#include <typeinfo>
#include <chrono>

#define F(x, y, z) (((x) & (y)) | ((~(x)) & (z))) //логические выражения для md5
#define G(x, y, z) (((x) & (z)) | ((y) & (~(z))))
#define H(x, y, z) ((x) ^ (y) ^ (z))
#define I(x, y, z) ((y) ^ ((x) | (~(z))))
#define LEFTROTATE(x, c) (((x) << (c)) | ((x) >> (32 - (c)))) // Левый циклический сдвиг для MD5
#define LOWBORDER  4
#define HIGHBORDER 1000 //константы
#define ZERO 0
#define MAXJACCARD 1
#define HASHSIZE 16
using namespace std;

//Класс для обработки и вывода исключений пользователю
class Exception : public exception
{
protected:
    char* str;
public:
    Exception(const char* s)
    {
        try {
            str = new char[strlen(s) + 1];
            strcpy_s(str, strlen(s) + 1, s);
        }
        catch (bad_alloc& e)
        {
            // если не удалось выделить память для str, то выводим сообщение об ошибке
            cout << "Memory allocation failed for str; " << e.what() << endl;
        }
    }
    Exception(const Exception& e)
    {
        try {
            str = new char[strlen(e.str) + 1];
            strcpy_s(str, strlen(e.str) + 1, e.str);
        }
        catch (bad_alloc& e)
        {
            // если не удалось выделить память для str, то выводим сообщение об ошибке
            cout << "Memory allocation failed for str; " << e.what() << endl;
        }
    }
    ~Exception()
    {
        delete[] str;
    }
    virtual void print() const
    {
        cout << "Exception: " << str << "; " << what();
    }
};

class WrongJaccardCoef : public Exception
{
public:
    WrongJaccardCoef(const char* s) : Exception(s) {}
    virtual void print() const
    {
        cout << "WrongJaccardCoef; " << what();
    }
};

class WrongHash : public Exception
{
public:
    WrongHash(const char* s) : Exception(s) {}
    virtual void print() const
    {
        cout << "WrongHash; " << what();
    }
};

template <typename T>
// Функция, которая вычисляет коэффициент Жаккара для двух векторов
double jaccard_coefficient(const vector<T>& v1, const vector<T>& v2) {
    // Находим пересечение векторов
    vector<T> intersection;
    set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(intersection));
    // Находим объединение векторов
    vector<T> union_;
    set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(union_));
    // Возвращаем отношение размера пересечения к размеру объединения
    double res = (double)intersection.size() / union_.size();
    if (res >= ZERO && res <= MAXJACCARD) return res;
    else throw WrongJaccardCoef("\nJaccar coeff less than 0 or more than 1");
}

// Функция MD5
string md5(string& message) {

    constexpr uint32_t s[64] = {
        7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7,
        12, 17, 22, 5, 9 , 14, 20, 5, 9 , 14, 20, 5,
        9 , 14, 20, 5, 9 , 14, 20, 4, 11, 16, 23, 4,
        11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23,
        6 ,10 ,15 ,21 ,6 ,10 ,15 ,21 ,6 ,10 ,15 ,21 ,
        6 ,10 ,15 ,21
    }; // Константы сдвига, определяющие количество бит, на которое нужно циклически сдвинуть каждое слово в четырех раундах алгоритма

    constexpr uint32_t K[64] = {
        0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee,
        0xf57c0faf, 0x4787c62a, 0xa8304613, 0xfd469501,
        0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be,
        0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821,
        0xf61e2562, 0xc040b340, 0x265e5a51, 0xe9b6c7aa,
        0xd62f105d, 0x02441453, 0xd8a1e681, 0xe7d3fbc8,
        0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed,
        0xa9e3e905, 0xfcefa3f8, 0x676f02d9, 0x8d2a4c8a,
        0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c,
        0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70,
        0x289b7ec6, 0xeaa127fa, 0xd4ef3085, 0x04881d05,
        0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665,
        0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039,
        0x655b59c3, 0x8f0ccc92, 0xffeff47d, 0x85845dd1,
        0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1,
        0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391
    }; // Константы смешивания, полученные из синуса целых чисел


    // Инициализация переменных
    uint32_t a0 = 0x67452301U; // A начальное значение первого регистра

    uint32_t b0 = 0xefcdab89U; // B начальное значение второго регистра

    uint32_t c0 = 0x98badcfeU; // C начальное значение третьего регистра

    uint32_t d0 = 0x10325476U; // D начальное значение четвертого регистра


    // Добавление бита "1" к сообщению
    message += char(0x80);

    // Добавление битов "0" до длины сообщения mod512 ==448
    while ((message.size() * 8) % 512 != 448) message += char(0); // Добавляем байты 00000000


    // Добавление длины сообщения в битах в конец сообщения
    uint64_t messageSizeBits = message.size() * 8; // Вычисляем длину сообщения в битах
    message.resize(message.size() + 8); // Увеличиваем размер сообщения на 8 байтов
    memcpy(&message[message.size() - 8], &messageSizeBits, 8);

    // Обработка каждого блока по 512 бит (64 байта)
    for (size_t i = 0; i < message.size(); i += 64)
    {
        array<uint32_t, 16> M; // Массив слов M[0..15]
        for (int j = 0; j < 16; j++) {
            M[j] = (uint8_t(message[i + j * 4]) << 0) // Собираем слово из четырех байтов в порядке от младшего к старшему
                | (uint8_t(message[i + j * 4 + 1]) << 8)
                | (uint8_t(message[i + j * 4 + 2]) << 16)
                | (uint8_t(message[i + j * 4 + 3]) << 24);
        }
        // Инициализация регистров для данного блока
        uint32_t A = a0;
        uint32_t B = b0;
        uint32_t C = c0;
        uint32_t D = d0;

        for (int j = 0; j < 64; j++)
        {
            uint32_t F, g;

            switch (j / 16)
            {
            case 0:
                F = F(B, C, D); // Функция F для первого раунда
                g = j; // Индекс слова M для первого раунда
                break;
            case 1:
                F = G(B, C, D); // Функция G для второго раунда
                g = (5 * j + 1) % 16; // Индекс слова M для второго раунда
                break;
            case 2:
                F = H(B, C, D); // Функция H для третьего раунда
                g = (3 * j + 5) % 16; // Индекс слова M для третьего раунда
                break;
            case 3:
                F = I(B, C, D); // Функция I для четвертого раунда
                g = (7 * j) % 16; // Индекс слова M для четвертого раунда
                break;
            }

            F += A + K[j] + M[g]; // Добавляем к F значение A, константу K и слово M
            A = D; // Перемещаем значение D в A
            D = C; // Перемещаем значение C в D
            C = B; // Перемещаем значение B в C
            B += LEFTROTATE(F, s[j]); // Прибавляем к B циклически сдвинутое значение F на s[j] бит влево
        }

        a0 += A; // Прибавляем к a0 значение A после 64 итераций
        b0 += B; // Прибавляем к b0 значение B после 64 итераций
        c0 += C; // Прибавляем к c0 значение C после 64 итераций
        d0 += D; // Прибавляем к d0 значение D после 64 итераций
    }

    // Формирование хеша из четырех 32-битных слов
    string hash;
    hash += char((a0 >> (0 * 8)) & 0xFF); // Добавляем к хешу младший байт из a0
    hash += char((a0 >> (1 * 8)) & 0xFF);
    hash += char((a0 >> (2 * 8)) & 0xFF);
    hash += char((a0 >> (3 * 8)) & 0xFF);

    hash += char((b0 >> (0 * 8)) & 0xFF); // Добавляем к хешу младший байт из b0
    hash += char((b0 >> (1 * 8)) & 0xFF);
    hash += char((b0 >> (2 * 8)) & 0xFF);
    hash += char((b0 >> (3 * 8)) & 0xFF);

    hash += char((c0 >> (0 * 8)) & 0xFF); // Добавляем к хешу младший байт из c0
    hash += char((c0 >> (1 * 8)) & 0xFF);
    hash += char((c0 >> (2 * 8)) & 0xFF);
    hash += char((c0 >> (3 * 8)) & 0xFF);

    hash += char((d0 >> (0 * 8)) & 0xFF); // Добавляем к хешу младший байт из d0
    hash += char((d0 >> (1 * 8)) & 0xFF);
    hash += char((d0 >> (2 * 8)) & 0xFF);
    hash += char((d0 >> (3 * 8)) & 0xFF);
    if (!hash.empty() && hash.size() == HASHSIZE)
        return hash;
    else
        throw WrongHash("\nEmpty or invalid hash");

}
template<typename T>
// Функция, которая генерирует хеш-функцию с помощью MD5 и соли
function<int(T)> md5_hash_function(const string& salt)
{
    return [&salt](const T& x) {
        string data;
        // if constexpr - стандарт C++ 17
        if constexpr (is_arithmetic_v<T>) {
            data = to_string(x);
        }
        else data = x;
        // Добавляем соль к данным
        data += salt;
        // Вызываем функцию MD5
        string hash = md5(data);
        // Преобразуем 128-битный хеш в 32-битное целое число
        uint32_t result = 0;
        for (int i = 0; i < 4; ++i)
        {
            result ^= (uint8_t(hash[i * 4]) << 0) //Применим XOR к каждому байту хеша
                | (uint8_t(hash[i * 4 + 1]) << 8) //Сдвигаем байты на соотвтствующее количество байт
                | (uint8_t(hash[i * 4 + 2]) << 16)
                | (uint8_t(hash[i * 4 + 3]) << 24);
        }
        if (result != ZERO) return result;
        else throw WrongHash("\nEmpty or invalid hash");
    };
}

template<typename T>
// Функция, которая вычисляет minhash для вектора данных с помощью k хеш-функций
vector<int> minhash(const vector<T>& data, int k) {
    vector<int> result(k, INT_MAX); // Результат - вектор из k минимальных хешей, инициализированный максимальным целым числом
    for (int i = 0; i < k; i++) {
        string salt = to_string(i);
        auto h = md5_hash_function<T>(salt); // Генерируем хеш-функцию с помощью MD5 и соли, равной номеру хеш-функции
        for (const T& x : data) {
            int hash = h(x); // Вычисляем хеш
            result[i] = result[i] ^ ((hash ^ result[i]) & -(hash < result[i]));// min через побитовые операции
        }
    }
    if (!result.empty()) return result;
    else throw WrongHash("\nEmpty or invalid hash");
}

// Функция, которая оценивает коэффициент Жаккара для двух векторов данных с помощью minhash
template<typename T>
double minhash_jaccard(const vector<T>& v1, const vector<T>& v2, int k) {
    // Вычисляем minhash для обоих векторов
    if (!v1.empty() && !v2.empty() && k > LOWBORDER && k < HIGHBORDER)
    {
        vector<int> m1 = minhash<T>(v1, k);
        vector<int> m2 = minhash<T>(v2, k);
        // Считаем, сколько раз минимальные хеши совпадают
        int matches = 0;
        for (int i = 0; i < k - 3; i += 4)
        {
            matches += (m1[i] == m2[i]) ? 1 : 0;
            matches += (m1[i + 1] == m2[i + 1]) ? 1 : 0;
            matches += (m1[i + 2] == m2[i + 2]) ? 1 : 0;
            matches += (m1[i + 3] == m2[i + 3]) ? 1 : 0;

        }
        return (double)matches / k; // Возвращаем отношение совпадений к количеству хеш-функций
    }

    else throw WrongHash("\nEmpty files or invalid count of hash function");
}


// Функция, которая выводит вектор на экран
template<typename T>
void print_vector(const vector<T>& v) {
    cout << "[ ";
    for (const T& i : v) cout << i << ' ';
    cout << "]\n";
}
// Вывод из файла в вектор
template<class T>
istream& operator>>(istream& s, vector<T>& v)
{
    if (typeid(s) == typeid(ifstream)) {
        T n;
        while (s >> n) {
            v.push_back(n);
        }
    }
    return s;
}
//Функция заполнения вектора данными из файла
template<class T>
void input(vector<T>& v, const string& filename) {
    ifstream file(filename);
    if (file) file >> v;
    else throw Exception("\nFile is wrong");
    file.close();
}

int main() {
    try {
        // Примеры векторов данных
        vector<int> v1 = { 1, 2, 3, 4, 5 };
        vector<int> v2 = { 2, 4, 6, 8, 10 };
        vector<int> v3 = { 1, 3, 5, 7, 9 };
        vector<int> v4 = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

        int k = 200; // Количество хеш-функций
        // Тестируем код на примерах
        cout << "Vector 1: ";
        print_vector(v1);
        cout << "Vector 2: ";
        print_vector(v2);
        cout << "Vector 3: ";
        print_vector(v3);
        cout << "Vector 4: ";
        print_vector(v4);


        cout << "Coef Jaccard for vectors 1 & 2: " << jaccard_coefficient(v1, v2) << "\n";
        cout << "Coef Jaccard for vectors 1 & 3: " << jaccard_coefficient(v1, v3) << "\n";
        cout << "Coef Jaccard for vectors 2 & 3: " << jaccard_coefficient(v2, v3) << "\n";
        cout << "Coef Jaccard for vectors 1 & 4: " << jaccard_coefficient(v1, v4) << "\n";
        cout << "Coef Jaccard for vectors 2 & 4: " << jaccard_coefficient(v2, v4) << "\n";


        cout << "Mark for coeff Jaccard with minhash for vectors 1&2: " << minhash_jaccard(v1, v2, k) << "\n";
        cout << "Mark for coeff Jaccard with minhash for vectors 1&3: " << minhash_jaccard(v1, v3, k) << "\n";
        cout << "Mark for coeff Jaccard with minhash for vectors 2&3: " << minhash_jaccard(v2, v3, k) << "\n";
        cout << "Mark for coeff Jaccard with minhash for vectors 1&4: " << minhash_jaccard(v1, v4, k) << "\n";
        cout << "Mark for coeff Jaccard with minhash for vectors 2&4: " << minhash_jaccard(v2, v4, k) << "\n";

        cout << "\n\n\n";
        vector<string> vs1; input(vs1, "big_text1.txt");
        vector<string> vs2; input(vs2, "big_text2.txt");
        cout << "Mark for coeff Jaccard with minhash for vectors(txt file) 1&2: " << minhash_jaccard(vs1, vs2, k) << "\n\n";
    }
    catch (Exception e)
    {
        cout << "\nCaught exeption: ", e.print();
        cout << "\n";
    }
    return 0;
}
