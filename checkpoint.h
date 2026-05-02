#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

// Định nghĩa ký tự phân cách khối (Block Separator)
const char CKP_SEPARATOR = '!';

class Checkpoint : public std::map<std::string, std::string> {
private:
    std::string checkpoint_filename;
    std::string current_prefix; // Lưu prefix hiện tại (VD: "StopRule!")

    Checkpoint() : current_prefix("") {}
    ~Checkpoint() {}
    Checkpoint(const Checkpoint&) = delete;
    Checkpoint& operator=(const Checkpoint&) = delete;

    std::string trim(const std::string& str);

    // Xử lý key: Nếu đang trong block, tự động nối thêm prefix
    std::string getFullKey(const std::string& key);

public:
    static Checkpoint& getInstance() {
        static Checkpoint instance;
        return instance;
    }

    void setFilename(const std::string& filename);

    // --- XỬ LÝ BLOCK ---
    void startBlock(const std::string& name);
    void endBlock();

    // --- CORE IO ---
    bool read(); 
    void dump(); 

    // --- CÁC HÀM GET ---
    int getInt(const std::string& key, int default_val = 0);
    double getDouble(const std::string& key, double default_val = 0.0);
    bool getBool(const std::string& key, bool default_val = false);
    std::string getString(const std::string& key, const std::string& default_val = "");

    // --- CÁC HÀM PUT ---
    void putInt(const std::string& key, int val);
    void putDouble(const std::string& key, double val);
    void putBool(const std::string& key, bool val);
    void putString(const std::string& key, const std::string& val);
};

#endif // CHECKPOINT_H