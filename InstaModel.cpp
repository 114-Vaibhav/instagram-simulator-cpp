#include <bits/stdc++.h>
#include <chrono>

using namespace std;

int noOfUsers = 500;
vector<float> eigenvector(noOfUsers, 0);
vector<float> katzCent(noOfUsers, 0);
vector<float> pagerank(noOfUsers, 0);
unordered_map<string, int> mp; // names to indices
vector<string> usernames;
vector<vector<int>> followMatrix(noOfUsers, vector<int>(noOfUsers, 0));
vector<vector<int>> postByUsers(noOfUsers);
vector<int> userByPost;
vector<string> posts;
unordered_map<string, int> postIndexMap;
vector<set<int>> likesOnPost;
vector<set<pair<string, int>>> commentOnPost;
vector<int> community(noOfUsers);
vector<set<int>> communities(10);

// read names from 'names.csv'
void readNamesFromFile(const string &filename)
{
    ifstream nameFile(filename);
    string line;
    int idx = 0;
    getline(nameFile, line); // Skip header line
    while (getline(nameFile, line))
    {
        usernames.push_back(line); // Store username
        mp[line] = idx++;          // Store username and corresponding index
    }
    nameFile.close();
}

void addANewCommunity()
{
    communities.push_back({});
}

void userJoinsCommunity()
{
    string you = "";
    cout << "Enter you username: ";
    cin >> you;
    while (mp.find(you) == mp.end())
    {
        cout << "Wrong username! Enter again: ";
        cin >> you;
    }
    int node = mp[you];
    int num;
    cout << "Enter from 0 to " << communities.size() - 1 << " to join a community: ";
    cin >> num;
    while (num < 0 || num >= communities.size())
    {
        cout << "Invalid community! Enter again: ";
        cin >> num;
    }
    communities[num].insert(node);
}

// read from 'edges.csv'
void readEdgesFromFile(const string &filename)
{
    ifstream edgeFile(filename);
    if (!edgeFile.is_open())
    {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }

    string line, user1, user2;
    getline(edgeFile, line); // Skip header line

    while (getline(edgeFile, line))
    {
        stringstream ss(line);
        getline(ss, user1, ',');
        getline(ss, user2);

        int u1_index = mp[user1];
        int u2_index = mp[user2];

        followMatrix[u1_index][u2_index] = 1; // Mark user1 follows user2
    }

    edgeFile.close();
}

// read from 'posts.csv'
void readPostsFromFile(const string &filename)
{
    ifstream postFile(filename);
    string line, post_name, user_name;
    int post_id = 0;         // Unique post ID (index)
    getline(postFile, line); // Skip header line
    while (getline(postFile, line))
    {
        stringstream ss(line);
        getline(ss, post_name, ',');
        getline(ss, user_name);
        posts.push_back(post_name);                 // Store post name
        postIndexMap[post_name] = post_id;          // Map post_name -> post_index
        int user_index = mp[user_name];             // Get user index from 'mp'
        userByPost.push_back(user_index);           // Store the user who posted this post
        postByUsers[user_index].push_back(post_id); // Map the post to the user (user index -> post index)
        post_id++;                                  // Increment post index for the next post
    }
    postFile.close();
}

void assignRandomCommunities(int numUsers)
{
    srand(time(0)); // Seed the random number generator

    // Assign a random community value (between 1 and 10)
    for (int i = 0; i < numUsers; i++)
    {
        community[i] = rand() % 10; // Random community index (0 to 9)
        communities[community[i]].insert(i);
    }
}

// read from 'likes.csv'
void readLikesFromFile(const string &filename)
{
    ifstream likeFile(filename);

    // Check if the file opened successfully
    if (!likeFile.is_open())
    {
        cout << "Error: Could not open likes.csv\n";
        return;
    }
    else
    {
        cout << "Successfully opened likes.csv\n";
    }

    // Resize likesOnPost to match the number of posts
    likesOnPost.resize(posts.size());

    string line, post_name, liker_name;

    // Read the header line and ensure it’s skipped
    if (!getline(likeFile, line))
    {
        cout << "Error: Empty file or failed to read header\n";
        return;
    }

    cout << "Starting to read likes data...\n";

    // Reading each line in the likes file
    while (getline(likeFile, line))
    {

        stringstream ss(line);
        getline(ss, post_name, ',');
        getline(ss, liker_name);

        // Retrieve indices and add the liker to the post’s like set
        int post_index = postIndexMap[post_name];
        int liker_index = mp[liker_name];

        likesOnPost[post_index].insert(liker_index);

        // Confirm that post_index is within bounds after resizing
        if (post_index >= likesOnPost.size())
        {
            cout << "Error: Post index " << post_index << " out of bounds for likesOnPost vector.\n";
            continue;
        }
    }

    likeFile.close();
}

//  Read from 'comments.csv'
void readCommentsFromFile(const string &filename)
{
    commentOnPost.resize(posts.size());
    ifstream commentFile(filename);
    string line, commenter_name, post_name, comment;
    getline(commentFile, line); // Skip header line
    while (getline(commentFile, line))
    {
        stringstream ss(line);
        getline(ss, post_name, ',');
        getline(ss, commenter_name, ',');
        getline(ss, comment);
        int post_index = postIndexMap[post_name]; // Get post index from postIndexMap
        int commenter_index = mp[commenter_name]; // Get commenter index

        commentOnPost[post_index].insert({comment, commenter_index}); // Add comment to the post's comment set
    }
    commentFile.close();
}

void dfs(vector<vector<vector<int>>> &kGraphs, vector<vector<int>> &graph, int begin, int node, int dist, int maxLevel)
{
    if (dist > maxLevel)
        return;

    if (dist >= 1)
        kGraphs[dist - 1][begin][node] = 1;

    for (int i = 0; i < graph[node].size(); i++)
    {
        if (graph[node][i] == 1)
        {
            dfs(kGraphs, graph, begin, i, dist + 1, maxLevel);
        }
    }
}

vector<float> katzCentrality(vector<vector<int>> &graph, vector<string> &usernames)
{
    int v = graph.size();
    vector<float> kc(v, 0);

    int k;
    cout << "Enter the max level: ";
    cin >> k;

    float alpha;
    cout << "Enter the damping factor : ";
    cin >> alpha;

    vector<vector<vector<int>>> kGraphs(k, vector<vector<int>>(v, vector<int>(v, 0)));

    for (int i = 0; i < v; i++)
    {
        dfs(kGraphs, graph, i, i, 0, k);
    }

    float ans;
    for (int l = 0; l < v; l++)
    {
        ans = 0;
        for (int i = 1; i <= k; i++)
        {
            for (int j = 0; j < v; j++)
            {
                ans += pow(alpha, i) * kGraphs[i - 1][j][l];
            }
        }
        kc[l] = ans;
    }

    for (int i = 0; i < kc.size(); i++)
    {
        cout << usernames[i] << " : " << kc[i] << endl;
    }

    return kc;
}

vector<int> djikstra(int from, vector<vector<int>> &graph)
{
    int v = graph.size();
    vector<int> distance(v, INT_MAX);
    distance[from] = 0;
    queue<pair<int, int>> st;
    st.push({0, from});
    while (st.size())
    {
        auto node = st.front();
        st.pop();

        int num = node.second, dist = node.first;
        for (auto it : graph[num])
        {
            if (distance[it] > distance[num] + 1)
            {
                st.push({distance[num] + 1, it});
                distance[it] = distance[num] + 1;
            }
        }
    }
    return distance;
}

void addANewAccount(vector<vector<int>> &followMatrix, unordered_map<string, int> &mp, int &numberOfUsers, vector<string> &usernames, vector<vector<int>> &postByUsers, vector<float> &eigenvector)
{
    string name;
    int n = followMatrix.size();
    cout << "Enter the new account username: ";
    cin >> name;
    while (mp.find(name) != mp.end())
    {
        cout << "Username already exists! Enter new username: ";
        cin >> name;
    }
    mp[name] = n;
    numberOfUsers++;
    usernames.push_back(name);
    followMatrix.push_back({});
    for (int i = 0; i <= n; i++)
    {
        if (i != n)
            followMatrix[i].push_back(0);
        followMatrix[n].push_back(0);
    }
    postByUsers.push_back({});
    eigenvector.push_back(0);
    katzCent.push_back(0);
    pagerank.push_back(0);
}

void AUnfollowsB(unordered_map<string, int> &mp, vector<vector<int>> &followMatrix)
{
    string A, B;
    int aIndex, bIndex;
    cout << "Enter the A's account's username: ";
    cin >> A;
    while (mp.find(A) == mp.end())
    {
        cout << "Username doesn't exist! Enter again: ";
        cin >> A;
    }
    aIndex = mp[A];
    cout << "Enter the account's username which get unfollow: ";
    cin >> B;
    while (mp.find(B) == mp.end())
    {
        cout << "Username doesn't exist! Enter again: ";
        cin >> B;
    }
    bIndex = mp[B];

    followMatrix[aIndex][bIndex] = 0;
}

void AfollowsB(unordered_map<string, int> &mp, vector<vector<int>> &followMatrix)
{
    string A, B;
    int aIndex, bIndex;
    cout << "Enter the account's username which followed: ";
    cin >> A;
    while (mp.find(A) == mp.end())
    {
        cout << "Username doesn't exist! Enter again: ";
        cin >> A;
    }
    aIndex = mp[A];
    cout << "Enter the account's username which got followed: ";
    cin >> B;
    while (mp.find(B) == mp.end())
    {
        cout << "Username doesn't exist! Enter again: ";
        cin >> B;
    }
    bIndex = mp[B];

    followMatrix[aIndex][bIndex] = 1;
}

void AveragePathLength(vector<vector<int>> &graph)
{
    int v = graph.size();
    vector<vector<int>> distances(v, vector<int>(v, INT_MAX));
    vector<vector<int>> AdjacencyList(v);
    for (int i = 0; i < v; i++)
    {
        distances[i][i] = 0;
        for (int j = 0; j < v; j++)
        {
            if (graph[i][j])
                AdjacencyList[i].emplace_back(j);
        }
    }

    for (int i = 0; i < v; i++)
    {
        distances[i] = djikstra(i, AdjacencyList);
    }
    long long total_dist = 0;
    for (int i = 0; i < v; i++)
    {
        for (int j = 0; j < v; j++)
        {
            if (distances[i][j] == INT_MAX)
            {
                cout << "All users are not connected!";
                return;
            }
            total_dist += distances[i][j];
        }
    }
    cout << "Average Path length is: " << (float)(total_dist) / (v * v) << endl;
    return;
}

void densityOfGraph(vector<vector<int>> &graph)
{
    int nodes = 0;
    int v = graph.size();
    for (int i = 0; i < v; i++)
    {
        for (int j = 0; j < v; j++)
        {
            if (graph[i][j])
                nodes++;
        }
    }
    cout << "The density of Network is: " << (float)(nodes) / (v * (v - 1)) << endl;
    return;
}

void reciprocity(vector<vector<int>> &graph)
{
    int v = graph.size();
    int biConnect = 0;
    for (int i = 0; i < v; i++)
    {
        for (int j = i + 1; j < v; j++)
        {
            if (graph[i][j] && graph[j][i])
                biConnect++;
        }
    }
    cout << "The number of reciprocal nodes are : " << biConnect;
}

void friends(unordered_map<string, int> &mp, vector<vector<int>> &graph, vector<string> &usernames)
{
    int v = graph.size();
    string username = "";
    cout << "Enter the username whose friends you want to find out: " << endl;
    cin >> username;
    while (mp.find(username) == mp.end())
    {
        cout << "Wrong username! Enter again: ";
        cin >> username;
    }
    int node = mp[username];
    cout << "Friends of " << username << " are:" << endl;
    for (int j = 0; j < v; j++)
    {
        if (graph[node][j] && graph[j][node])
            cout << usernames[j] << endl;
    }
}

void addLikesOnPost(vector<set<int>> &likesOnPosts, unordered_map<string, int> &mp)
{
    string you = "", whose = "";
    cout << "Enter you username: ";
    cin >> you;
    while (mp.find(you) == mp.end())
    {
        cout << "Wrong username! Enter again: ";
        cin >> you;
    }
    int node = mp[you];

    cout << "Enter whose post you want to like: ";
    cin >> whose;
    while (mp.find(whose) == mp.end())
    {
        cout << "Wrong username! Enter again: ";
        cin >> whose;
    }
    int node2 = mp[whose];

    int postId;
    for (auto it : postByUsers[node2])
    {
        cout << it << " " << posts[it] << endl;
    }
    cout << "Enter the post id you want to like: ";
    cin >> postId;
    while (postId >= likesOnPosts.size())
    {
        cout << "Post with this id doesn't exist! Enter again: ";
        cin >> postId;
    }
    likesOnPosts[postId].insert(node);
}

void addCommentOnPost(vector<set<pair<string, int>>> &commentOnPosts, unordered_map<string, int> &mp)
{
    string you = "", whose = "";
    cout << "Enter you username: ";
    cin >> you;
    while (mp.find(you) == mp.end())
    {
        cout << "Wrong username! Enter again: ";
        cin >> you;
    }
    int node = mp[you];

    cout << "Enter whose post you want to comment on: ";
    cin >> whose;
    while (mp.find(whose) == mp.end())
    {
        cout << "Wrong username! Enter again: ";
        cin >> whose;
    }
    int node2 = mp[whose];
    int postId;
    for (auto it : postByUsers[node2])
    {
        cout << it << " : " << posts[it] << endl;
    }

    cout << "Enter the post id you want to comment on : ";
    cin >> postId;
    while (postId >= commentOnPosts.size())
    {
        cout << "Post with this id doesn't exist! Enter again: ";
        cin >> postId;
    }
    string comment;
    cout << "Enter your comment: ";
    cin.ignore(1, '/n');
    getline(cin, comment);
    commentOnPosts[postId].insert({comment, node});
}

void addNewPost(vector<int> &userByPost, vector<vector<int>> &postByUsers, unordered_map<string, int> &mp, vector<set<int>> &likesOnPosts, vector<set<pair<string, int>>> &commentOnPosts)
{
    int newPostNo = userByPost.size();
    string uploader = "";
    cout << "Enter username of who want to upload post: ";
    cin >> uploader;
    while (mp.find(uploader) == mp.end())
    {
        cout << "Wrong username! Enter correct username: ";
        cin >> uploader;
    }
    int node = mp[uploader];
    string postMessage = "";
    cout << "Enter the post message: ";
    cin.ignore(1, '/n');
    getline(cin, postMessage);
    userByPost.push_back(node);
    posts.push_back(postMessage);
    postByUsers[node].push_back(newPostNo);
    likesOnPosts.push_back({});
    commentOnPosts.push_back({});
}

void whoLiked(vector<set<int>> &likesOnPosts, vector<string> &usernames)
{
    string you = "";
    cout << "Enter username whose post you want to see: ";
    cin >> you;
    while (mp.find(you) == mp.end())
    {
        cout << "Wrong username! Enter again: ";
        cin >> you;
    }
    int node = mp[you];

    for (auto it : postByUsers[node])
    {
        cout << it << " : " << posts[it] << endl;
    }
    int post;
    cout << "Enter the post id :" << endl;
    cin >> post;
    while (post >= likesOnPosts.size())
    {
        cout << "Wrong post id!";
        cin >> post;
    }
    cout << "The people who liked post " << post << " are:";
    for (auto it : likesOnPosts[post])
    {
        cout << usernames[it] << ",";
    }
}

void whoCommented(vector<set<pair<string, int>>> &commentOnPost, vector<string> &usernames)
{
    string you = "";
    cout << "Enter username on whose post you want to see comments: ";
    cin >> you;
    while (mp.find(you) == mp.end())
    {
        cout << "Wrong username! Enter again: ";
        cin >> you;
    }
    int node = mp[you];

    for (auto it : postByUsers[node])
    {
        cout << it << " : " << posts[it] << endl;
    }
    int post;
    cout << "Enter the post id :" << endl;
    cin >> post;
    while (post >= commentOnPost.size())
    {
        cout << "Wrong post id!" << endl;
        cin >> post;
    }
    cout << "The comments on post " << post << " are: " << endl;
    for (auto it : commentOnPost[post])
    {
        cout << "\"" << it.first << "\" : " << usernames[it.second] << endl;
    }
}

void LCC(vector<vector<int>> &graph, unordered_map<string, int> &mp)
{
    int v = graph.size();
    string user = "";
    cout << "Enter the username to find his/her LCC: ";
    cin >> user;
    while (mp.find(user) == mp.end())
    {
        cout << "Wrong username! Enter correct username: ";
        cin >> user;
    }
    int check_node = mp[user];
    vector<int> neighbours;
    for (int i = 0; i < v; i++)
    {
        if (graph[check_node][i])
        {
            neighbours.push_back(i);
        }
    }
    int total = 0, common = 0;
    for (int i = 0; i < neighbours.size(); i++)
    {
        for (int j = i + 1; j < neighbours.size(); j++)
        {
            int x = neighbours[i], y = neighbours[j];
            if (graph[x][y])
                common++;
            total++;
        }
    }
    cout << "The LCC  of " << user << " is:" << (float)common / total;
}

void GCC(vector<vector<int>> &graph)
{
    int v = graph.size(), triplets = 0, total = 0;
    for (int i = 0; i < v; i++)
    {
        for (int j = i + 1; j < v; j++)
        {
            for (int k = j + 1; k < v; k++)
            {
                if (graph[i][j] && graph[j][k] && graph[k][i])
                    triplets++;
                total++;
            }
        }
    }
    cout << "The GCC of model is: " << (float)triplets / total;
    return;
}

float localClusteringCoefficient(int check_node, vector<vector<int>> &graph)
{
    int v = graph.size();
    vector<int> neighbours;
    for (int i = 0; i < v; i++)
    {
        if (graph[check_node][i] || graph[i][check_node])
        {
            neighbours.push_back(i);
        }
    }
    int total = 0, common = 0;
    for (int i = 0; i < neighbours.size(); i++)
    {
        for (int j = 0; j < neighbours.size(); j++)
        {
            if (i == j)
                continue;
            int x = neighbours[i], y = neighbours[j];
            if (graph[x][y] || graph[y][x])
                common++;
            total++;
        }
    }
    return (float)common / total;
}

void clusterCoeffIfTriadic(vector<vector<int>> &graph, unordered_map<string, int> &mp)
{
    int v = graph.size();
    int node1, node2;
    string user1, user2;

    cout << "Enter first user: ";
    cin >> user1;
    while (mp.find(user1) == mp.end())
    {
        cout << "Wrong username! Enter correct username: ";
        cin >> user1;
    }
    cout << "Enter second user: ";
    cin >> user2;
    while (mp.find(user2) == mp.end())
    {
        cout << "Wrong username! Enter correct username: ";
        cin >> user2;
    }
    node1 = mp[user1];
    node2 = mp[user2];

    float ans = 0;
    for (int i = 0; i < v; i++)
    {
        if (graph[i][node1] && graph[i][node2])
        {
            float temp = localClusteringCoefficient(i, graph);
            ans = max(ans, temp);
        }
    }
    cout << "Probability of node between " << user1 << " and " << user2 << " is " << ans;
}

// Matrix Multiplication

vector<float> matrixMultiply(vector<float> eigen, vector<vector<int>> &graph)
{

    int v = eigen.size();
    vector<float> values(v);
    float normalizedValue = 0;
    for (int i = 0; i < v; i++)
    {
        float sum = 0;
        for (int j = 0; j < v; j++)
        {
            sum += eigen[j] * graph[j][i];
        }
        normalizedValue += pow(sum, 2);
        values[i] = sum;
    }
    normalizedValue = pow(normalizedValue, 0.5);
    for (int i = 0; i < v; i++)
    {
        values[i] /= normalizedValue;
    }
    return values;
}

vector<float> eigenVector(vector<vector<int>> graph, vector<string> &usernames)
{
    int v = graph.size();
    for (int i = 0; i < v; i++)
    {
        for (int j = i + 1; j < v; j++)
        {
            swap(graph[i][j], graph[j][i]);
        }
    }
    int iterate;
    cout << "Enter the number of iterations you want: ";
    cin >> iterate;

    vector<float> eigen(v, 1);
    for (int i = 0; i < iterate; i++)
    {
        eigen = matrixMultiply(eigen, graph);
        cout << endl;
    }
    for (int i = 0; i < v; i++)
    {
        cout << usernames[i] << " : " << eigen[i] << endl;
    }
    return eigen;
}

// Converts edge list to Adjacency Matrix
void convertToMatrix(vector<vector<int>> &edges, vector<vector<int>> &followMatrix)
{
    int n = edges.size();
    for (int i = 0; i < n; i++)
    {
        int x = edges[i][0], y = edges[i][1];
        followMatrix[x][y] = 1;
    }
}

// to add usernames to map and assign them id
void inserts(vector<string> &usernames, unordered_map<string, int> &mp)
{
    for (int i = 0; i < usernames.size(); i++)
    {
        mp[usernames[i]] = i;
    }
}

void printMat(vector<vector<float>> &mat, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << fixed << setprecision(3) << mat[i][j] << " ";
        }
        cout << endl;
    }
}

vector<vector<float>> NormAdj(vector<vector<int>> &adj, int n)
{
    vector<vector<float>> norm(n, vector<float>(n, 0));
    for (int i = 0; i < n; i++)
    {
        int cnt = 0;
        for (int j = 0; j < n; j++)
        {
            if (adj[i][j] == 1)
                cnt++;
        }
        for (int j = 0; j < n; j++)
        {
            if (adj[i][j] == 1)
                norm[i][j] = float(1.0 / cnt);
        }
    }
    return norm;
}

void print1D(vector<float> &x, int n, vector<string> &userNames)
{
    cout << endl;
    for (int i = 0; i < n; i++)
    {
        cout << userNames[i] << " " << fixed << setprecision(3) << x[i] << endl;
    }
    cout << endl;
}

// multiply 2 matrices
vector<float> mult(vector<float> &R0, vector<vector<float>> &norm, int n)
{
    vector<float> res(n);
    for (int j = 0; j < n; j++)
    {
        float ans = 0.0;
        for (int i = 0; i < n; i++)
        {
            ans += 0.8 * norm[i][j] * R0[i];
        }
        res[j] = ans;
    }
    return res;
}

// add two vectors
vector<float> add(vector<float> &a, vector<float> &b, int n)
{
    vector<float> ans(n);
    for (int i = 0; i < n; i++)
    {
        ans[i] = a[i] + b[i];
    }
    return ans;
}

vector<float> pageRank(vector<vector<int>> &adjL, int n, vector<string> &userNames)
{
    vector<vector<float>> norm = NormAdj(adjL, n);
    float d = 0.8;
    vector<float> E(n, (1 - d) * (1.0 / n));

    vector<float> R0(n, 1.0 / n);

    vector<float> secondPart_R1 = mult(R0, norm, n);
    vector<float> pageRank_R1 = add(E, secondPart_R1, n);
    cout << "\n......Printing Page Rank R1.....\n";
    print1D(pageRank_R1, n, userNames);

    vector<float> secondPart_R2 = mult(pageRank_R1, norm, n);
    vector<float> pageRank_R2 = add(E, secondPart_R2, n);
    cout << "\n......Printing Page Rank R2.....\n";
    print1D(pageRank_R2, n, userNames);

    vector<float> secondPart_R3 = mult(pageRank_R2, norm, n);
    vector<float> pageRank_R3 = add(E, secondPart_R3, n);
    cout << "\n......Printing Page Rank R3.....\n";
    print1D(pageRank_R3, n, userNames);

    return pageRank_R3;
}

void celebrityFollowing(vector<vector<int>> &followMatrix, vector<string> &usernames, unordered_map<string, int> &mp)
{
    string name = "";
    cout << "Enter your username: ";
    cin >> name;
    while (mp.find(name) == mp.end())
    {
        cout << "Username does not exist in Network! Enter again: ";
        cin >> name;
    }
    int node = mp[name];
    int n = followMatrix.size();
    vector<int> followers(n, 0);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (followMatrix[j][i])
                followers[i]++;
        }
    }
    priority_queue<pair<int, int>> pq;
    for (int i = 0; i < followers.size(); i++)
    {
        pq.push({followers[i], i});
    }
    cout << "The celebrities " << name << " follows are: ";
    for (int i = 0; i < 5 && pq.size(); i++)
    {
        auto temps = pq.top();
        pq.pop();
        cout << usernames[temps.second] << " : " << temps.first << endl;
    }
}

void friendSuggestion(vector<vector<int>> &followMatrix, vector<string> &usernames, unordered_map<string, int> &mp)
{
    priority_queue<pair<int, string>> pq;
    string name = "";
    int n = followMatrix.size();
    cout << "Enter the username for which you want friend suggestions: ";
    cin >> name;
    while (mp.find(name) == mp.end())
    {
        cout << "Wrong username! Enter again :";
        cin >> name;
    }
    int node = mp[name];
    for (int i = 0; i < n; i++)
    {
        if (i == node)
            continue;
        int sameFollowings = 0;
        for (int j = 0; j < n; j++)
        {
            if (followMatrix[node][j] && followMatrix[i][j])
            {
                sameFollowings++;
            }
        }
        if (sameFollowings)
            pq.push({sameFollowings, usernames[i]});
    }
    cout << "The top recommendations for " << name << " are: " << endl;
    for (int i = 0; i < 5 && pq.size(); i++)
    {
        auto nums = pq.top();
        pq.pop();
        cout << nums.second << " -> " << nums.first << " common followings" << endl;
    }
}

vector<pair<int, int>> findShortPaths(int ptr, vector<vector<int>> &AM)
{
    int v = AM.size();
    vector<pair<int, int>> Paths(v, {INT_MAX, 0});
    // first stores distance, second stores no of paths

    queue<pair<int, int>> q;
    q.push({ptr, -1});
    int dist = 0;
    while (q.size())
    {
        int x = q.size();

        for (int i = 0; i < x; i++)
        {
            int node = q.front().first;
            int par = q.front().second;
            // cout<<node<<" "<<par<<Paths[node].first<<dist<<endl;
            q.pop();
            if (Paths[node].first < dist)
                continue;
            else if (Paths[node].first == dist)
            {
                if (par == -1)
                    Paths[node] = {0, 1};
                else
                    Paths[node] = {dist, Paths[node].second + Paths[par].second};
            }
            else
            {
                Paths[node].first = dist;
                if (par != -1)
                    Paths[node].second = Paths[par].second;
                else
                    Paths[node].second = 1;

                for (int j = 0; j < v; j++)
                {
                    if (AM[node][j] == 1)
                        q.push({j, node});
                }
            }
        }
        dist++;
    }

    return Paths;
}

void displayData()
{
    // Display usernames and their indices
    cout << "Usernames and their indices:\n";
    for (const auto &entry : mp)
    {
        cout << "User: " << entry.first << ", Index: " << entry.second << endl;
    }

    // Display the follow matrix
    cout << "\nFollow Matrix (Adjacency Matrix):\n";
    for (int i = 0; i < noOfUsers; i++)
    {
        for (int j = 0; j < noOfUsers; j++)
        {
            cout << followMatrix[i][j] << " ";
        }
        cout << endl;
    }

    // Display posts and the user who posted them
    cout << "\nPosts and their users:\n";
    for (size_t i = 0; i < posts.size(); i++)
    {
        cout << "Post: " << posts[i] << ", Posted by: " << usernames[userByPost[i]] << endl;
    }

    // Display likes on each post
    cout << "\nLikes on posts:\n";
    for (size_t i = 0; i < likesOnPost.size(); i++)
    {
        cout << "Post: " << posts[i] << " - Liked by users: ";
        for (int user : likesOnPost[i])
        {
            cout << usernames[user] << " ";
        }
        cout << endl;
    }

    // Display comments on each post
    cout << "\nComments on posts:\n";
    for (size_t i = 0; i < commentOnPost.size(); i++)
    {
        cout << "Post: " << posts[i] << " - Comments: ";
        for (const auto &comment : commentOnPost[i])
        {
            cout << "(" << comment.first << ", " << usernames[comment.second] << ") ";
        }
        cout << endl;
    }
}

float ratioGraph(int node, vector<vector<pair<int, int>>> &shortestPaths)
{
    int v = shortestPaths.size();
    float ans = 0;
    for (int i = 0; i < v; i++)
    {
        if (i == node || shortestPaths[i][node].first == INT_MAX)
            continue;
        for (int j = 0; j < v; j++)
        {
            if (j == node || shortestPaths[node][j].first == INT_MAX)
                continue;
            if (shortestPaths[i][j].first == shortestPaths[i][node].first + shortestPaths[node][j].first)
            {
                float temp = shortestPaths[i][node].second * shortestPaths[node][j].second;
                temp = temp / shortestPaths[i][j].second;
                ans = ans + temp;
            }
        }
    }

    return ans;
}

void betweenCent(vector<string> &usernames, vector<vector<int>> &AM)
{
    cout << endl;
    int v = AM.size();
    vector<vector<pair<int, int>>> shortestPaths(v, vector<pair<int, int>>(v, {INT_MAX, 0}));
    for (int i = 0; i < v; i++)
    {
        shortestPaths[i] = findShortPaths(i, AM);
    }

    vector<float> btwnCentrality(v, -1);

    for (int i = 0; i < v; i++)
    {
        btwnCentrality[i] = ratioGraph(i, shortestPaths);
        cout << usernames[i] << " : " << btwnCentrality[i] << endl;
    }
}

void searchSuggestion()
{
    string name = "";
    cout << "Enter the username to search: ";
    cin >> name;
    vector<int> similars;
    for (int i = 0; i < usernames.size(); i++)
    {
        if (usernames[i].find(name) != -1)
        {
            similars.push_back(i);
        }
    }
    if (similars.size() == 0)
    {
        cout << "No similar username!";
        return;
    }
    sort(similars.begin(), similars.end(), [&](int i, int j)
         {
        if(pagerank[i]>=pagerank[j]) return true;
        return false; });
    cout << "The most oprimal searches for this name are: " << endl;
    for (int i = 0; i < similars.size(); i++)
    {
        cout << usernames[similars[i]] << endl;
    }
}

void importantUserByPostWeightage()
{
    priority_queue<pair<int, int>> postWeight;

    for (int i = 0; i < noOfUsers; i++)
    {
        int weight = 0;
        for (int j = 0; j < postByUsers[i].size(); j++)
        {
            weight = weight + likesOnPost[postByUsers[i][j]].size() + 2 * commentOnPost[postByUsers[i][j]].size();
        }
        if (postByUsers[i].size() == 0)
        {
            postWeight.push({0, i});
            continue;
        }
        weight /= postByUsers[i].size();
        postWeight.push({weight, i});
    }
    int i = 0;
    cout << "The top accounts based on engagement are: " << endl;
    while (i < 5 && postWeight.size())
    {
        auto details = postWeight.top();
        postWeight.pop();
        int userNode = details.second;
        string username = usernames[userNode];
        cout << username << endl;
        i++;
    }
}

int calculateIndegreeWithinCommunity(int userID, const vector<int> &community, const vector<vector<int>> &followMatrix)
{
    int indegree = 0;
    int communityID = community[userID];

    for (int i = 0; i < followMatrix.size(); i++)
    {
        if (followMatrix[i][userID] == 1 && community[i] == communityID)
        {
            indegree++;
        }
    }

    return indegree;
}

int calculateOutdegree(int userID, const vector<int> &community, const vector<vector<int>> &followMatrix)
{
    int outdegree = 0;
    int communityID = community[userID];

    for (int i = 0; i < followMatrix.size(); i++)
    {
        if (followMatrix[userID][i] == 1 && community[i] != communityID)
        {
            outdegree++;
        }
    }

    return outdegree;
}

double calculateInternalClusteringCoefficient(int userID, const vector<int> &community, const vector<vector<int>> &followMatrix)
{
    int communityID = community[userID];
    set<int> neighborsInCommunity;

    for (int i = 0; i < followMatrix.size(); i++)
    {
        if (followMatrix[userID][i] == 1 && community[i] == communityID)
        {
            neighborsInCommunity.insert(i);
        }
    }

    int actualEdges = 0;
    for (int neighbor : neighborsInCommunity)
    {
        for (int otherNeighbor : neighborsInCommunity)
        {
            if (neighbor != otherNeighbor && followMatrix[neighbor][otherNeighbor] == 1)
            {
                actualEdges++;
            }
        }
    }

    actualEdges /= 2;

    int possibleEdges = neighborsInCommunity.size() * (neighborsInCommunity.size() - 1) / 2;

    if (possibleEdges == 0)
    {
        return 0.0;
    }

    return static_cast<double>(actualEdges) / possibleEdges;
}

double calculatePermanence(int userID, const vector<int> &community, const vector<vector<int>> &followMatrix)
{
    int indegree = calculateIndegreeWithinCommunity(userID, community, followMatrix);
    int degree = 0;

    for (int i = 0; i < followMatrix.size(); i++)
    {
        if (followMatrix[userID][i] == 1 || followMatrix[i][userID] == 1)
        {
            degree++;
        }
    }

    int outdegree = calculateOutdegree(userID, community, followMatrix);
    double clusteringCoefficient = calculateInternalClusteringCoefficient(userID, community, followMatrix);

    double permanence = static_cast<double>(indegree) / (degree * outdegree) - (1 - clusteringCoefficient);

    return permanence;
}

double calculateModularity(const vector<vector<int>> &followMatrix, const vector<int> &community)
{
    int numUsers = followMatrix.size();
    int numCommunities = 10;

    double m = 0;
    for (int i = 0; i < numUsers; ++i)
    {
        for (int j = i + 1; j < numUsers; ++j)
        {
            if (followMatrix[i][j] == 1)
            {
                m += 1;
            }
        }
    }

    double modularity = 0;

    for (int i = 0; i < numUsers; ++i)
    {
        for (int j = 0; j < numUsers; ++j)
        {
            if (followMatrix[i][j] == 1)
            {
                int communityI = community[i];
                int communityJ = community[j];

                int degreeI = 0, degreeJ = 0;
                for (int k = 0; k < numUsers; ++k)
                {
                    if (followMatrix[i][k] == 1)
                        degreeI++;
                    if (followMatrix[j][k] == 1)
                        degreeJ++;
                }

                double contribution = (1.0 / (2 * m)) * (followMatrix[i][j] - (degreeI * degreeJ) / (2 * m));

                if (communityI == communityJ)
                {
                    modularity += contribution;
                }
            }
        }
    }

    return modularity;
}

int calculateDegree(int userID, const vector<vector<int>> &followMatrix)
{
    int degree = 0;
    int numUsers = followMatrix.size();
    for (int i = 0; i < numUsers; ++i)
    {
        if (followMatrix[userID][i] == 1)
        {
            degree++;
        }
    }
    return degree;
}

double calculateGeneralizedPermanence(int userID, const vector<vector<int>> &followMatrix,
                                      const vector<int> &community)
{
    int numUsers = followMatrix.size();

    int indegreeWithinCommunity = 0;
    int degreeOfNode = calculateDegree(userID, followMatrix);

    for (int i = 0; i < numUsers; ++i)
    {
        if (followMatrix[i][userID] == 1 && community[i] == community[userID])
        {
            indegreeWithinCommunity++;
        }
    }

    int connectionsOutsideCommunity = 0;
    for (int i = 0; i < numUsers; ++i)
    {
        if (followMatrix[i][userID] == 1 && community[i] != community[userID])
        {
            connectionsOutsideCommunity++;
        }
    }

    double internalClusteringCoefficient = calculateInternalClusteringCoefficient(userID, community, followMatrix);

    double permanence = (double)indegreeWithinCommunity / (connectionsOutsideCommunity * degreeOfNode) - (1 - internalClusteringCoefficient);

    return permanence;
}

void printRelevantData(vector<string> &usernames, vector<vector<int>> &followMatrix, vector<vector<int>> &postByUsers)
{
    int v = usernames.size();
    for (int i = 0; i < v; i++)
    {

        int followers = 0, followings = 0;
        for (int j = 0; j < v; j++)
        {
            if (followMatrix[i][j])
                followings++;
            if (followMatrix[j][i])
                followers++;
        }
        cout << "Name: " << usernames[i] << endl;
        cout << "No. of following accounts: " << followings << endl;
        cout << "No. of followers: " << followers << endl;
        cout << "No. of posts: " << postByUsers[i].size() << endl;
        cout << "--------------------------------------------------------" << endl;
    }
}

int main()
{

    readNamesFromFile("names.csv");
    noOfUsers = usernames.size();
    cout << "hellp" << endl;
    noOfUsers = usernames.size();
    readEdgesFromFile("edges80K.csv");
    cout << "hello" << endl;

    readPostsFromFile("posts.csv");
    cout << "Hell" << endl;
    readLikesFromFile("likes.csv");
    cout << "hee;" << endl;
    readCommentsFromFile("comments.csv");

    assignRandomCommunities(noOfUsers);

    int option = -2;
    while (option != -1)
    {
        cout << endl;
        cout << "Enter your choice:" << endl;
        cout << "-1 to EXIT" << endl;
        cout << "0. To show info about everyone: " << endl;
        // Time Complexity: O(n)
        cout << "1. To add a new user to model: " << endl;
        // Time Complexity: O(n)
        cout << "2. User A starts following user B" << endl;
        // Time Complexity: O(1)
        cout << "3. To calculate Average Path Length: " << endl;
        // Time Complexity: O(n^3)
        cout << "4. Density of graph : " << endl;
        // Time Complexity: O(n^2)
        cout << "5. Number of reciprocal Nodes: " << endl;
        // Time Complexity: O(n^2)
        cout << "6. Find No. of friends in groups: " << endl;
        // Time Complexity: O(n)
        cout << "7. Check LCC of a user: " << endl;
        // Time Complexity: O(n^3)
        cout << "8. Check GCC of the network: " << endl;
        // Time Complexity: O(n^2)
        cout << "9. Check if an edge will exist between 2 users: " << endl;
        // Time Complexity: O(n^3)
        cout << "10. Find eigenvector value : " << endl;
        // Time Complexity: O(n^3)
        cout << "11. Find Katz centrality of nodes: " << endl;
        // Time Complexity: O(n^3)
        cout << "12. Find PageRank : " << endl;
        // Time Complexity: O(n^2*k)
        cout << "13. Add a new Post: " << endl;
        // Time Complexity: O(1)
        cout << "14. Like a post : " << endl;
        // Time Complexity: O(1)
        cout << "15. Comment on a post:" << endl;
        // Time Complexity: O(1)
        cout << "16. Who all liked a post: " << endl;
        // Time Complexity: O(n)
        cout << "17. Comments of a post : " << endl;
        // Time Complexity: O(n)
        cout << "18. How many celebrities a person follows: " << endl;
        // Time Complexity: O(n^2)
        cout << "19. For checking friend suggestions: " << endl;
        // Time Complexity: O(n^2)
        cout << "20. Find Betweenness Centrality: " << endl;
        // Time Complexity: O(n^3)
        cout << "21. User A unfollows user B: " << endl;
        // Time Complexity: O(1)
        cout << "22. Search if username not clear: " << endl;
        // Time Complexity: O(1)
        cout << "23. To find which users' posts have more engagement: " << endl;
        // Time Complexity: O(n)
        cout << "24. Find Permanence " << endl;
        // Time Complexity: O(n^2)
        cout << "25. Find Modularity " << endl;
        // Time Complexity: O(n^3)
        cout << "26. Find Gen Permanence " << endl;
        // Time Complexity: O(n^2)
        cout << "27. Add a new community " << endl;
        // Time Complexity: O(1)
        cout << "28. User joins a new community " << endl;

        cout << "Enter your choice: ";
        cin >> option;

        auto start = std::chrono::high_resolution_clock::now();

        if (option == 0)
            printRelevantData(usernames, followMatrix, postByUsers);
        else if (option == 1)
            addANewAccount(followMatrix, mp, noOfUsers, usernames, postByUsers, eigenvector);
        else if (option == 2)
            AfollowsB(mp, followMatrix);
        else if (option == 3)
            AveragePathLength(followMatrix);
        else if (option == 4)
            densityOfGraph(followMatrix);
        else if (option == 5)
            reciprocity(followMatrix);
        else if (option == 6)
            friends(mp, followMatrix, usernames);
        else if (option == 7)
            LCC(followMatrix, mp);
        else if (option == 8)
            GCC(followMatrix);
        else if (option == 9)
            clusterCoeffIfTriadic(followMatrix, mp);
        else if (option == 10)
            eigenvector = eigenVector(followMatrix, usernames);
        else if (option == 11)
            katzCent = katzCentrality(followMatrix, usernames);
        else if (option == 12)
            pagerank = pageRank(followMatrix, noOfUsers, usernames);
        else if (option == 13)
            addNewPost(userByPost, postByUsers, mp, likesOnPost, commentOnPost);
        else if (option == 14)
            addLikesOnPost(likesOnPost, mp);
        else if (option == 15)
            addCommentOnPost(commentOnPost, mp);
        else if (option == 16)
            whoLiked(likesOnPost, usernames);
        else if (option == 17)
            whoCommented(commentOnPost, usernames);
        else if (option == 18)
            celebrityFollowing(followMatrix, usernames, mp);
        else if (option == 19)
            friendSuggestion(followMatrix, usernames, mp);
        else if (option == 20)
            betweenCent(usernames, followMatrix);
        else if (option == 21)
            AUnfollowsB(mp, followMatrix);
        else if (option == 22)
            searchSuggestion();
        else if (option == 23)
            importantUserByPostWeightage();
        else if (option == 24)
        {

            string name;
            cout << "Enter User name: ";
            cin >> name;
            cout << " " << calculatePermanence(mp[name], community, followMatrix);
        }
        else if (option == 25)
        {
            cout << "Modularity is : " << calculateModularity(followMatrix, community) << endl;
        }
        else if (option == 26)
        {
            string name;
            cout << "Enter User name: ";
            cin >> name;

            cout << " " << calculateGeneralizedPermanence(mp[name], followMatrix, community);
        }
        else if (option == 27)
        {
            addANewCommunity();
        }
        else if (option == 28)
        {
            userJoinsCommunity();
        }

        else if (option == -1)
            break;
        else if (option != -1)
        {
            cout << "Invalid Option!";
            cin >> option;
        }

        auto stop = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

        cout << "Time taken: " << duration.count() << " ms" << endl;
    }
    return 0;
}
