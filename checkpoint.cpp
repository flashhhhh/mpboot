#include "checkpoint.h"
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cctype>

// ==========================================
// CÁC HÀM TIỆN ÍCH (UTILITIES)
// ==========================================
std::string Checkpoint::trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (std::string::npos == first) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, (last - first + 1));
}

std::string Checkpoint::getFullKey(const std::string& key) {
    // Nếu key do user truyền vào đã có sẵn dấu '!', không nối thêm prefix nữa
    if (key.find(CKP_SEPARATOR) != std::string::npos) {
        return key;
    }
    return current_prefix + key;
}

void Checkpoint::startBlock(const std::string& name) {
    current_prefix = name + CKP_SEPARATOR;
}

void Checkpoint::endBlock() {
    current_prefix = "";
}

void Checkpoint::setFilename(const std::string& filename) {
    checkpoint_filename = filename;
}

// ==========================================
// THUẬT TOÁN ĐỌC FILE (FLATTEN TO MAP)
// ==========================================
bool Checkpoint::read() {
    std::ifstream in(checkpoint_filename);
    if (!in.is_open()) return false;

    this->clear();
    std::string line;
    std::string file_prefix = ""; 

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line.find("---") == 0) continue;

        // KIỂM TRA THỤT LỀ: Nếu dòng không có khoảng trắng ở đầu -> Thoát khỏi block
        size_t indent = line.find_first_not_of(" \t");
        if (indent == 0) {
            file_prefix = ""; 
        }

        std::string trimmed_line = trim(line);
        size_t colon_pos = trimmed_line.find(':');
        
        if (colon_pos != std::string::npos) {
            std::string key = trim(trimmed_line.substr(0, colon_pos));
            std::string value = trim(trimmed_line.substr(colon_pos + 1));
            
            if (value.empty()) {
                // Nếu không có value, đây là Header của Block (VD: "StopRule:")
                file_prefix = key + CKP_SEPARATOR;
            } else {
                // Lưu vào RAM dưới dạng "StopRule!curIteration" = "102"
                (*this)[file_prefix + key] = value;

                printf("Key: %s, value: %s\n", (file_prefix + key).c_str(), value.c_str());
            }
        }
    }
    in.close();
    return true;
}

// ==========================================
// THUẬT TOÁN GHI FILE (UNFLATTEN FROM MAP)
// ==========================================
void Checkpoint::dump() {
    if (checkpoint_filename.empty()) return;

    std::string tmp_file = checkpoint_filename + ".tmp";
    std::ofstream out(tmp_file, std::ios::out | std::ios::trunc);
    if (!out.is_open()) return;

    out << "--- # MPBoot Checkpoint\n";

    std::string last_prefix = "";

    // Lưu ý: std::map sẽ tự động sort các key theo thứ tự Alphabet.
    // Điều này có nghĩa là "StopRule!curIteration" và "StopRule!time_vec" chắc chắn sẽ nằm cạnh nhau.
    for (const auto& pair : *this) {
        std::string full_key = pair.first;
        std::string value = pair.second;
        
        size_t sep_pos = full_key.find(CKP_SEPARATOR);
        
        if (sep_pos != std::string::npos) {
            // Đây là một biến nằm trong Block
            std::string prefix = full_key.substr(0, sep_pos);
            std::string local_key = full_key.substr(sep_pos + 1);
            
            // Nếu đây là lần đầu ghi biến của Block này, in tên Block ra trước
            if (prefix != last_prefix) {
                out << prefix << ":\n";
                last_prefix = prefix;
            }
            // In biến có thụt lề
            out << "  " << local_key << ": " << value << "\n";
        } else {
            // Đây là biến Global
            last_prefix = ""; 
            out << full_key << ": " << value << "\n";
        }
    }

    out.close();
    std::rename(tmp_file.c_str(), checkpoint_filename.c_str());
}

// ==========================================
// CÁC HÀM TRUY XUẤT GET / PUT ĐÃ HOÀN THIỆN
// ==========================================

// --- Integer ---
int Checkpoint::getInt(const std::string& key, int default_val) {
    auto it = this->find(getFullKey(key));
    if (it != this->end() && !it->second.empty()) {
        try { return std::stoi(it->second); } catch(...) {}
    }
    return default_val;
}

void Checkpoint::putInt(const std::string& key, int val) {
    (*this)[getFullKey(key)] = std::to_string(val);
}

// --- Double ---
double Checkpoint::getDouble(const std::string& key, double default_val) {
    auto it = this->find(getFullKey(key));
    if (it != this->end() && !it->second.empty()) {
        try { return std::stod(it->second); } catch(...) {}
    }
    return default_val;
}

void Checkpoint::putDouble(const std::string& key, double val) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(9) << val; // Giữ độ chính xác cao cho tính toán cây
    (*this)[getFullKey(key)] = oss.str();
}

// --- Boolean ---
bool Checkpoint::getBool(const std::string& key, bool default_val) {
    auto it = this->find(getFullKey(key));
    if (it != this->end() && !it->second.empty()) {
        std::string val = it->second;
        // Chuyển về chữ thường để so sánh cho an toàn (True, TRUE, true đều pass)
        std::transform(val.begin(), val.end(), val.begin(),
            [](unsigned char c){ return std::tolower(c); });
        return (val == "true" || val == "1");
    }
    return default_val;
}

void Checkpoint::putBool(const std::string& key, bool val) {
    (*this)[getFullKey(key)] = val ? "true" : "false";
}

// --- String ---
std::string Checkpoint::getString(const std::string& key, const std::string& default_val) {
    auto it = this->find(getFullKey(key));
    if (it != this->end()) {
        return it->second;
    }
    return default_val;
}

void Checkpoint::putString(const std::string& key, const std::string& val) {
    (*this)[getFullKey(key)] = val;
}