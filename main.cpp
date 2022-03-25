#include<bits/stdc++.h>
using namespace std;

using q_elem = pair<int, double>;

struct Compare{
    bool operator() (q_elem const a, q_elem const b) {
        if(a.second < b.second) {
            return true;
        }
        return (a.second == b.second) && a.first < b.first;
    }
};

using dataQueue = priority_queue<q_elem, vector<q_elem>, Compare>;

double cosine_dist(vector<double> a, vector<double> b){
    assert(a.size() == b.size());
    int n = a.size();
    double dotProd = 0.0;
    double lenA = 0.0;
    double lenB = 0.0;
    for(int i = 0; i < n; i++) {
        double ad = a[i];
        double bd = b[i];
        dotProd += (ad * bd);
        lenA += (ad * ad);
        lenB += (bd * bd);
    }
    dotProd = 1.0 - dotProd / (sqrt(lenA) * sqrt(lenB));
    return dotProd;
}

void trim(dataQueue &dq, int limit){
    while (dq.size() > limit)
        dq.pop();
}

dataQueue maxToMinOrViceVersaHeap(dataQueue dq) {
    vector<q_elem> v;
    dataQueue newDq;
    while (!dq.empty()) {
        q_elem top = dq.top();
        newDq.push({-top.first, -top.second});
        dq.pop();
    }
    return newDq;
}

void SearchLayer(vector<double> q, int k, vector<uint32_t> indptr, vector<int32_t> index, vector<uint32_t> level_offset, int currLevel, unordered_set<int> &visited, vector<vector<double>> vect, dataQueue &candidates) {
    dataQueue newCandidates = maxToMinOrViceVersaHeap(candidates);
    while(newCandidates.size() > 0) {
        int curr = -newCandidates.top().first;
        newCandidates.pop();
        int start = indptr[curr] + level_offset[currLevel];
        int end = indptr[curr] + level_offset[currLevel + 1];
        for(int i = start;i<end;i++){
            int currIdx = index[i];
            if(visited.count(currIdx) || currIdx == -1) {
                continue;
            }
            visited.insert(currIdx);
            double currDist = cosine_dist(q, vect[currIdx]);
            if (currDist > candidates.top().second)
                continue;
            candidates.push(make_pair(currIdx, currDist));
            trim(candidates, k);
            newCandidates.push(make_pair(-currIdx, -currDist));
        }
    }
}

dataQueue queryHNSW(vector<double> q, int top_k, int ep, vector<uint32_t> indptr, vector<int32_t> index, vector<uint32_t> level_offset, int max_level, vector<vector<double>> vect){
    dataQueue candidates;
    candidates.push(make_pair(ep, cosine_dist(q, vect[ep])));
    unordered_set<int> visited;
    visited.insert(ep);
    for(int i = max_level;i>=0;i--) {
        SearchLayer(q, top_k, indptr, index, level_offset, i, visited, vect, candidates);
    } 
    return candidates;
}

ifstream open_file(string filename){
    ifstream file;
    file.open(filename, ios::binary);
    if(!file.is_open()){
        cout << "Could not open file" << endl;
        exit(1);
    }
    return file;
}

template <typename T, size_t S>
void read_val(T *n, ifstream &file) {
    unsigned char buffer[S];
    file.read((char *)buffer, S);
    *n = *(T *)buffer;
}

int main(int argc, char *argv[]){
    assert(argc == 5);

    string data_dir = argv[1];
    int top_k = stoi(argv[2]);
    string user = argv[3];
    string out_file = argv[4];

    ifstream params_file = open_file(data_dir + "/params");
    uint32_t ep, max_level, l, d;
    read_val<uint32_t, 4>(&ep, params_file);
    read_val<uint32_t, 4>(&max_level, params_file);
    read_val<uint32_t, 4>(&l, params_file);
    read_val<uint32_t, 4>(&d, params_file);
    cout << "ep: " << ep << ' ' << "max_level: " << max_level << ' ' << "l: " << l << ' ' << "d: " << d << endl;

    ifstream indptr_file = open_file(data_dir + "/indptr");
    vector<uint32_t> indptr;
    while (indptr_file.peek() != EOF) {
        uint32_t curr;
        read_val<uint32_t, 4>(&curr, indptr_file);
        indptr.push_back(curr);
    }
    cout << "indptr read" << endl;

    ifstream index_file = open_file(data_dir + "/index");
    vector<int32_t> index;
    while (index_file.peek() != EOF) {
        int32_t curr;
        read_val<int32_t, 4>(&curr, index_file);
        index.push_back(curr);
    }
    cout << "index read" << endl;

    ifstream level_offset_file = open_file(data_dir + "/level_offset");
    vector<uint32_t> level_offset;
    while (level_offset_file.peek() != EOF) {
        uint32_t n;
        read_val<uint32_t, 4>(&n, level_offset_file);
        level_offset.push_back(n);
    }
    cout << "level_offset read" << endl;

    ifstream vect_file = open_file(data_dir + "/vect");
    vector<vector<double>> vect;
    while (vect_file.peek() != EOF) {
        vector<double> items;
        for (int i = 0; i < d; i++) {
            double dbl;
            read_val<double, 8>(&dbl, vect_file);
            items.push_back(dbl);
        }
        vect.push_back(items);
    }

    ifstream user_file(user);
    if (!user_file.is_open()) {
        cout << "Could not open file" << endl;
        exit(1);
    }
    vector<vector<double>> users;
    while (user_file.peek() != EOF) {
        vector<double> items;
        for (int i = 0; i < d; i++) {
            double dbl;
            user_file >> dbl;
            items.push_back(dbl);
        }
        users.push_back(items);
    }

    dataQueue candidates = queryHNSW(users[0], top_k, ep, indptr, index, level_offset, max_level, vect);

    while (!candidates.empty()) {
        cout << candidates.top().first << ' ' << candidates.top().second << endl;
        candidates.pop();
    }
}