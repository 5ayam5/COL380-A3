#include<bits/stdc++.h>
using namespace std;

struct Data{
    int numItems;
    vector<double> dataItems;

    Data(int num, vector<double> items){
        this->numItems = num;
        for(int i=0;i<num;i++){
            this->dataItems.push_back(items[num]);
        }
    } 
};

class Compare{
public:
    bool operator() (q_elem a, q_elem b){
        if(a.second <= b.second){
            return true;
        }else if(a.second == b.second){
            return a.first <= b.first;
        }else{
            return false;
        }
    }
};

typedef pair<int, double> q_elem;
typedef priority_queue<q_elem, vector<q_elem>, Compare> dataQueue;

double cosine_dist(Data a, Data b){
    assert(a.numItems == b.numItems);
    double dotProd = 0.0;
    double lenA = 0.0;
    double lenB = 0.0;
    for(int i=0;i<a.numItems;i++){
        double ad = a.dataItems[i];
        double bd = b.dataItems[i];
        dotProd += (ad * bd);
        lenA += (ad * ad);
        lenB += (bd * bd);
    }
    dotProd = 1.0 - dotProd / (sqrt(lenA) * sqrt(lenB));
    return dotProd;
}

void trim(dataQueue &dq, int limit){
    vector<q_elem> v;
    while(dq.size()>0){
        v.push_back(dq.top());
        dq.pop();
    }
    for(int i=0;i<limit;i++){
        dq.push(v[i]);
    }
}

void SearchLayer(Data q, vector<uint32_t> indptr, vector<int32_t> index, vector<uint32_t> level_offset, int currLevel, unordered_set<int> &visited, vector<Data> vect, dataQueue &candidates){
    dataQueue *newCandidates = new dataQueue(candidates);
    int k = candidates.size();
    while(newCandidates->size() > 0){
        int curr = newCandidates->top().first;
        newCandidates->pop();
        int start = indptr[curr] + level_offset[currLevel];
        int end = indptr[curr] + level_offset[currLevel + 1];
        for(int i = start;i<end;i++){
            int currIdx = index[i];
            if(visited.find(currIdx) == visited.end() || currIdx == -1){
                continue;
            }
            visited.insert(currIdx);
            double currDist = cosine_dist(q, vect[currIdx]);
            //Do distance check
            candidates.push(make_pair(currIdx, currDist));
            trim(candidates,k);
            newCandidates->push(make_pair(currIdx,currDist));
        }
    }
}

dataQueue queryHNSW(Data q, int top_k, int ep, vector<uint32_t> indptr, vector<int32_t> index, vector<uint32_t> level_offset, int max_level, vector<Data> vect){
    dataQueue candidates;
    candidates.push(make_pair(ep, cosine_dist(q, vect[ep])));
    unordered_set<int> visited;
    visited.insert(ep);
    for(int i = max_level - 1;i>=0;i--){
        SearchLayer(q, indptr, index, level_offset, i, visited, vect, candidates);
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
void read_val(T &n, ifstream &file) {
    memset(&n, 0, sizeof(n));
    unsigned char buffer[S];
    file.read((char *)buffer, S);
    for (size_t i = 0; i < S; i++) {
        n = (n << 8) | buffer[i];
    }
}

int main(int argc, char *argv[]){
    assert(argc == 5);

    string data_dir = argv[1];
    int top_k = stoi(argv[2]);
    string user = argv[3];
    string out_file = argv[4];

    ifstream params_file = open_file(data_dir + "/params");
    uint32_t ep, max_level, l, d;
    read_val<uint32_t, 4>(ep, params_file);
    read_val<uint32_t, 4>(max_level, params_file);
    read_val<uint32_t, 4>(l, params_file);
    read_val<uint32_t, 4>(d, params_file);

    ifstream indptr_file = open_file(data_dir + "/indptr");
    vector<uint32_t> indptr;
    for(int i=0;i<l+1;i++){
        uint32_t curr;
        read_val<uint32_t, 4>(curr, indptr_file);
        indptr.push_back(curr);
    }

    ifstream index_file = open_file(data_dir + "/index");
    vector<int32_t> index;
    while (index_file.peek() != EOF) {
        int32_t n;
        read_val<int32_t, 4>(n, index_file);
        index.push_back(n);
    }

    ifstream level_offset_file = open_file(data_dir + "/level_offset");
    vector<uint32_t> level_offset;
    while (level_offset_file.peek() != EOF) {
        uint32_t n;
        read_val<uint32_t, 4>(n, level_offset_file);
        level_offset.push_back(n);
    }

    ifstream vect_file = open_file(data_dir + "/vect");
    vector<Data> vect;
    while (vect_file.peek() != EOF) {
        vector<double> items;
        for (int i = 0; i < d; i++) {
            double dbl;
            read_val<double, 8>(dbl, vect_file);
            items.push_back(dbl);
        }
        vect.push_back(Data(d, items));
    }

    ifstream user_file(user);
    if (!user_file.is_open()) {
        cout << "Could not open file" << endl;
        exit(1);
    }
    vector<Data> users;
    while (user_file.peek() != EOF) {
        vector<double> items;
        for (int i = 0; i < d; i++) {
            double dbl;
            user_file >> dbl;
            items.push_back(dbl);
        }
        users.push_back(Data(d, items));
    }

    queryHNSW(users[0], top_k, ep, indptr, index, level_offset, max_level, vect);
}