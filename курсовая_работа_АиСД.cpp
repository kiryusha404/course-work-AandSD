#include <iostream>
#include <vector>
#include <queue>
#include <climits>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <numeric> 
#include <iomanip>
#include <string>


using namespace std;

// Функция реализует алгоритм Прима
void primMST(vector<vector<int>> graph) {
    int v = graph.size();

    // вектор для хранения родителя вершины
    vector<int> parent(v);

    // вектор хранит вес/стоимость 
    vector<int> key(v);

    // вектор для представления набора
    // вершин, включённых в MST
    vector<bool> vis(v);

    priority_queue<pair<int, int>,
        vector<pair<int, int>>,
        greater<pair<int, int>>> pq;

    // Инициализируем все элементы вектора key как INFINITE
    // а элементы вектора vis — как false
    for (int i = 0; i < v; i++) {
        key[i] = INT_MAX;
        vis[i] = false;
    }

    // Всегда включаем первую вершину в минимальное оставное дерево
    // Устанавливаем ключ равным 0, чтобы эта вершина
    // была выбрана в качестве первой.
    key[0] = 0;

    // Первая вершина всегда является корнем минимального оставного дерева
    parent[0] = -1;

    // Добавляем исходную вершину в минимальную кучу
    pq.push({ 0, 0 });

    while (!pq.empty()) {
        int node = pq.top().second;
        pq.pop();
        vis[node] = true;
        for (int i = 0; i < v; i++) {

            // Если вершина не была посещена
            // и вес ребра соседней вершины меньше, чем значение ключа
            // соседней вершины, то обновляем его.
            if (!vis[i] && graph[node][i] != 0
                && graph[node][i] < key[i]) {
                pq.push({ graph[node][i], i });
                key[i] = graph[node][i];
                parent[i] = node;
            }
        }
    }

    // Выводим рёбра и их веса 
    //cout << "Алгоритм Прима\n";
    //cout << "Ребро \tВес\n";
    //for (int i = 1; i < v; i++) {
    //    cout << parent[i] << " - " << i
    //        << " \t" << graph[i][parent[i]] << " \n";
    //}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Структура для представления ребра
struct Edge {
    int u, v, weight;

    // Для сортировки по весу
    bool operator<(const Edge& other) const {
        return weight < other.weight;
    }
};

// Класс для системы непересекающихся множеств (Union-Find)
class UnionFind {
private:
    vector<int> parent;
    vector<int> rank;

public:
    UnionFind(int n) {
        parent.resize(n);
        rank.resize(n, 0);
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
        }
    }

    // Найти корень множества, содержащего x (с оптимизацией сжатия пути)
    int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);
        }
        return parent[x];
    }

    // Объединить множества, содержащие x и y (с оптимизацией по рангу)
    void unionSets(int x, int y) {
        int rootX = find(x);
        int rootY = find(y);

        if (rootX != rootY) {
            if (rank[rootX] < rank[rootY]) {
                parent[rootX] = rootY;
            }
            else if (rank[rootX] > rank[rootY]) {
                parent[rootY] = rootX;
            }
            else {
                parent[rootY] = rootX;
                rank[rootX]++;
            }
        }
    }
};

// Функция реализует алгоритм Краскала
void kruskalMST(const vector<vector<int>>& graph) {
    int n = graph.size();
    vector<Edge> edges;

    // 1. Преобразуем матрицу смежности в список рёбер
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {  // j = i+1: избегаем дублирования и петли
            if (graph[i][j] != 0) {        // есть ребро
                edges.push_back({ i, j, graph[i][j] });
            }
        }
    }

    // Сортируем рёбра по весу (по возрастанию)
    sort(edges.begin(), edges.end());

    // Инициализируем Union-Find для n вершин
    UnionFind uf(n);

    // Строим минимальное оставное дерево
    vector<Edge> mst;  // рёбра 
    int totalWeight = 0;

    for (const Edge& e : edges) {
        int u = e.u, v = e.v;
        // Если вершины в разных множествах — ребро не создаёт цикл
        if (uf.find(u) != uf.find(v)) {
            uf.unionSets(u, v);
            mst.push_back(e);
            totalWeight += e.weight;
        }
    }
    
    
    // Выводим рёбра и их веса 
    //cout << "Алгоритм Краскала\n";
    //cout << "Ребро \tВес\n";
    //for (const Edge& e : mst) {
    //    cout << e.u << " - " << e.v
    //        << " \t" << e.weight << " \n";
    //}
    //
}


// Функция для генерации случайной симметричной матрицы
vector<vector<int>> generateRandomSymmetricMatrix(int n) {
    srand(static_cast<unsigned int>(time(0)));
    vector<vector<int>> matrix(n, vector<int>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            int val = rand();
            matrix[i][j] = val;
            matrix[j][i] = val;
        }
    }
    return matrix;
}

int main() {
    setlocale(0, "");
    srand(static_cast<unsigned int>(time(0)));

    // Фиксированные ширины колонок (в символах)
    const int COL1_WIDTH = 14;  // Размер матрицы
    const int COL2_WIDTH = 16;  // Алгоритм Прима
    const int COL3_WIDTH = 16;  // Алгоритм Краскала

    // Шапка таблицы
    cout << std::setw(COL1_WIDTH) << std::left << "Размер матрицы";
    cout << " | ";
    cout << std::setw(COL2_WIDTH) << std::left << "Алгоритм Прима";
    cout << " | ";
    cout << std::setw(COL3_WIDTH) << std::left << "Алгоритм Краскала";
    cout << "\n";

    // Разделительная линия
    cout << std::string(COL1_WIDTH, '-');
    cout << "-|-";
    cout << std::string(COL2_WIDTH, '-');
    cout << "-|-";
    cout << std::string(COL3_WIDTH, '-');
    cout << "\n";

    int sizes[] = { 3, 5, 7, 10, 15, 20, 30, 50, 80, 100 };

    for (int check = 0; check < 10; ++check) {
        int n = sizes[check];
        vector<chrono::microseconds> durations_prim;
        vector<chrono::microseconds> durations_kruskal;

        // 1000 замеров
        for (int i = 0; i < 100000; ++i) {
            auto matrix = generateRandomSymmetricMatrix(n);

            // Прим
            auto start = chrono::high_resolution_clock::now();
            primMST(matrix);
            auto end = chrono::high_resolution_clock::now();
            durations_prim.push_back(
                chrono::duration_cast<chrono::microseconds>(end - start)
            );

            // Краскал
            start = chrono::high_resolution_clock::now();
            kruskalMST(matrix);
            end = chrono::high_resolution_clock::now();
            durations_kruskal.push_back(
                chrono::duration_cast<chrono::microseconds>(end - start)
            );
        }

        // Средние значения
        auto total_prim = accumulate(durations_prim.begin(), durations_prim.end(), chrono::microseconds(0));
        auto avg_prim = total_prim / 100000;

        auto total_kruskal = accumulate(durations_kruskal.begin(), durations_kruskal.end(), chrono::microseconds(0));
        auto avg_kruskal = total_kruskal / 100000;

        // Вывод строки таблицы
        cout << std::setw(COL1_WIDTH) << std::left << (to_string(n) + "x" + to_string(n));
        cout << " | ";
        cout << std::setw(COL2_WIDTH - 4) << std::right << avg_prim.count() << " мкс";  // -4: учитываем " мкс"
        cout << " | ";
        cout << std::setw(COL3_WIDTH - 4) << std::right << avg_kruskal.count() << " мкс";
        cout << "\n";
    }

    return 0;
}

