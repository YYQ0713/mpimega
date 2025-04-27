#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <map>
#include <sstream>
#include <algorithm>
#include <fcntl.h>
#include <sys/wait.h>
#include <thread>
#include <future>

#include "mpienv/mpienv.hpp"
#include "definitions.h"
#include "utils/cpu_dispatch.h"
#include "utils/utils.h"

int main_assemble(int argc, char **argv, MPIEnviroment &mpienv);
int main_local(int argc, char **argv);
int main_iterate(int argc, char **argv);
int main_build_lib(int argc, char **argv, MPIEnviroment &mpienv);

int main_kmer_count(int argc, char **argv, MPIEnviroment &mpienv);
int main_read2sdbg(int argc, char **argv, MPIEnviroment &mpienv);
int main_seq2sdbg(int argc, char **argv, MPIEnviroment &mpienv);

int main_contig2fastg(int argc, char **argv);
int main_read_stat(int argc, char **argv);
int main_filter_by_len(int argc, char **argv);

namespace fs = std::filesystem;

class Usage : public std::runtime_error {
public:
    explicit Usage(const std::string& message) : std::runtime_error(message) {}
};

std::string abspath(const std::string &path) {
    // 展开 "~"
    std::string expanded = path;
    if (!expanded.empty() && expanded[0] == '~') {
        if (const char* home = std::getenv("HOME")) {
            expanded = std::string(home) + expanded.substr(1);
        }
    }
    // 将路径转换为绝对路径，并进行词法归一化（去除冗余的"./"和"../"）
    fs::path p(expanded);
    fs::path absPath = fs::absolute(p).lexically_normal();
    return absPath.string();
}

void mkdir_if_not_exists(const std::string& file_name) {
    if (!fs::exists(file_name)) {
        try {
            fs::create_directories(file_name);
        } catch (const fs::filesystem_error& e) {
            std::cerr << "Error creating directory: " << e.what() << std::endl;
            throw Usage("Failed to create directory: " + file_name);
        }
    }
}


const std::string usage_message = R"(
contact: Dinghua Li <voutcn@gmail.com>

Usage:
  megahit [options] {-1 <pe1> -2 <pe2> | --12 <pe12> | -r <se>} [-o <out_dir>]

  Input options that can be specified for multiple times (supporting plain text and gz/bz2 extensions)
    -1                       <pe1>          comma-separated list of fasta/q paired-end #1 files, paired with files in <pe2>
    -2                       <pe2>          comma-separated list of fasta/q paired-end #2 files, paired with files in <pe1>
    --12                     <pe12>         comma-separated list of interleaved fasta/q paired-end files
    -r/--read                <se>           comma-separated list of fasta/q single-end files

Optional Arguments:
  Basic assembly options:
    --min-count              <int>          minimum multiplicity for filtering (k_min+1)-mers [2]
    --k-list                 <int,int,..>   comma-separated list of kmer size
                                            all must be odd, increment <= 28
                                            [21,29,39,59,79,99,119,141]

  Another way to set --k-list (overrides --k-list if one of them set):
    --k-min                  <int>          minimum kmer size, must be odd number [21]
    --k-max                  <int>          maximum kmer size, must be odd number [141]
    --k-step                 <int>          increment of kmer size of each iteration (<= 28), must be even number [10]

  Advanced assembly options:
    --no-mercy                              do not add mercy kmers
    --bubble-level           <int>          intensity of bubble merging (0-2), 0 to disable [2]
    --merge-level            <l,s>          merge complex bubbles of length <= l*kmer_size and similarity >= s [20,0.95]
    --prune-level            <int>          strength of low depth pruning (0-3) [2]
    --prune-depth            <int>          remove unitigs with avg kmer depth less than this value [2]
    --disconnect-ratio       <float>        disconnect unitigs if its depth is less than this ratio times 
                                            the total depth of itself and its siblings [0.1]  
    --low-local-ratio        <float>        remove unitigs if its depth is less than this ratio times
                                            the average depth of the neighborhoods [0.2]
    --max-tip-len            <int>          remove tips less than this value [2*k]
    --cleaning-rounds        <int>          number of rounds for graph cleanning [5]
    --no-local                              disable local assembly
    --kmin-1pass                            use 1pass mode to build SdBG of k_min

  Presets parameters:
    --presets                <str>          override a group of parameters; possible values:
                                            meta-sensitive: '--min-count 1 --k-list 21,29,39,49,...,129,141'
                                            meta-large: '--k-min 27 --k-max 127 --k-step 10'
                                            (large & complex metagenomes, like soil)

  Hardware options:
    -m/--memory              <float>        max memory in byte to be used in SdBG construction
                                            (if set between 0-1, fraction of the machine's total memory) [0.9]
    --mem-flag               <int>          SdBG builder memory mode. 0: minimum; 1: moderate;
                                            others: use all memory specified by '-m/--memory' [1]
    -t/--num-cpu-threads     <int>          number of CPU threads [# of logical processors]
    --no-hw-accel                           run MEGAHIT without BMI2 and POPCNT hardware instructions
    -n/--num-processes       <int>          mpirun with num of processes[1]

  Output options:
    -o/--out-dir             <string>       output directory [./megahit_out]
    --out-prefix             <string>       output prefix (the contig file will be OUT_DIR/OUT_PREFIX.contigs.fa)
    --min-contig-len         <int>          minimum length of contigs to output [200]
    --keep-tmp-files                        keep all temporary files
    --tmp-dir                <string>       set temp directory

Other Arguments:
    --continue                              continue a MEGAHIT run from its last available check point.
                                            please set the output directory correctly when using this option.
    --test                                  run MEGAHIT on a toy test dataset
    -h/--help                               print the usage message
    -v/--version                            print version
)";


std::string inpipe_cmd(const std::string& fileName) {
    if (fileName.size() > 3 && fileName.compare(fileName.size() - 3, 3, ".gz") == 0) {
        return "gzip -cd " + fileName;
    } else if (fileName.size() > 4 && fileName.compare(fileName.size() - 4, 4, ".bz2") == 0) {
        return "bzip2 -cd " + fileName;
    } else {
        return "";
    }
}

int execlp_cmd(const std::string& fileName) {
    if (fileName.size() > 3 && fileName.compare(fileName.size() - 3, 3, ".gz") == 0) {
        return 3;
    } else if (fileName.size() > 4 && fileName.compare(fileName.size() - 4, 4, ".bz2") == 0) {
        return 5;
    } else {
        return 0;
    }
}

// 检测系统可用内存，单位：字节
size_t detect_available_mem() {
    // 尝试使用 sysconf 获取页大小和物理页数
    long psize = sysconf(_SC_PAGE_SIZE);
    long pcount = sysconf(_SC_PHYS_PAGES);
    if (psize > 0 && pcount > 0) {
        return static_cast<size_t>(psize) * static_cast<size_t>(pcount);
    } else {
        // 如果 sysconf 无法获取，则根据平台采用不同方案
    #if defined(__APPLE__)
        // Mac OS X: 使用 sysctl 命令获取 hw.memsize
        FILE* pipe = popen("sysctl hw.memsize", "r");
        if (!pipe) {
            throw std::runtime_error("popen failed");
        }
        char buffer[128];
        std::string result;
        while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            result += buffer;
        }
        pclose(pipe);
        // 例如返回 "hw.memsize: 17179869184\n"
        std::istringstream iss(result);
        std::string label;
        size_t memsize = 0;
        if (!(iss >> label >> memsize)) {
            throw std::runtime_error("Failed to parse hw.memsize");
        }
        return memsize;
    #elif defined(__linux__)
        // Linux: 使用 free 命令
        FILE* pipe = popen("free", "r");
        if (!pipe) {
            throw std::runtime_error("popen failed");
        }
        char buffer[256];
        size_t mem_total = 0;
        // 跳过第一行标题
        if (fgets(buffer, sizeof(buffer), pipe) == nullptr) {
            pclose(pipe);
            throw std::runtime_error("Failed to read free output");
        }
        // 第二行通常形如 "Mem:  16342356  ..."
        if (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            std::istringstream iss(buffer);
            std::string key;
            if (!(iss >> key >> mem_total)) {
                pclose(pipe);
                throw std::runtime_error("Failed to parse free output");
            }
            // free 返回的是 KB，转换为字节
            mem_total *= 1024;
        }
        pclose(pipe);
        if (mem_total == 0) {
            throw std::runtime_error("Failed to detect available memory from free command");
        }
        return mem_total;
    #else
        throw std::runtime_error("Platform not supported for memory detection");
    #endif
    }
}

// 根据 kmer_k 构造图前缀，并保证相应目录存在
std::string graph_prefix(const std::string &temp_dir, int kmer_k) {
    // 拼接目录 temp_dir/k<kmer_k>
    fs::path k_dir = fs::path(temp_dir) / ("k" + std::to_string(kmer_k));
    mkdir_if_not_exists(k_dir);
    // 返回 temp_dir/k<kmer_k>/<kmer_k>
    fs::path prefix = k_dir / std::to_string(kmer_k);
    return prefix.string();
}

std::string contig_prefix(const std::string &contig_dir, int kmer_k) {
    // 拼接目录 contig_dir/k<kmer_k>
    fs::path prefix = fs::path(contig_dir) / ("k" + std::to_string(kmer_k));
    return prefix.string();
}

class Options {
public:
    std::string out_dir = "";
    std::string temp_dir = "";
    bool test_mode = false;
    bool continue_mode = false;
    bool force_overwrite = false;
    float memory = 0.9;
    int min_contig_len = 200;
    int k_min = 21;
    int k_max = 141;
    int k_step = 10;
    std::vector<int> k_list = {21, 29, 39, 59, 79, 99, 119, 141};
    bool auto_k = true;
    bool set_list_by_min_max_step = false;
    int min_count = 2;
    bool has_popcnt = true;
    bool hw_accel = true;
    int max_tip_len = -1;
    bool no_mercy = false;
    bool no_local = false;
    int bubble_level = 2;
    int merge_len = 20;
    float merge_similar = 0.95;
    int prune_level = 2;
    int prune_depth = 2;
    int num_cpu_threads = 0;
    int num_processes = 1;
    float disconnect_ratio = 0.1;
    float low_local_ratio = 0.2;
    int cleaning_rounds = 5;
    bool keep_tmp_files = false;
    int mem_flag = 1;
    std::string out_prefix = "";
    bool kmin_1pass = false;
    std::vector<std::string> pe1;
    std::vector<std::string> pe2;
    std::vector<std::string> pe12;
    std::vector<std::string> se;
    std::string presets = "";
    bool verbose = false;
    MPIEnviroment mpienv_;

    std::string log_file_name() {
        if (out_prefix.empty()) {
            return out_dir + "/log";
        } else {
            return out_dir + "/" + out_prefix + ".log";
        }
    }

    std::string option_file_name() {
        return out_dir + "/options.json";
    }

    std::string contig_dir() {
        return out_dir + "/intermediate_contigs";
    }

    std::string read_lib_path() {
        return temp_dir + "/reads.lib";
    }

    std::string get_executable_path() {

        char buffer[1024];
        ssize_t length = readlink("/proc/self/exe", buffer, sizeof(buffer) - 1);
        if (length == -1) {
            // 处理错误
            return "";
        }
        buffer[length] = '\0';
        return std::string(buffer);
    }

    size_t host_mem() {
        if (memory <= 0) {
            throw Usage("Please specify a positive number for -m flag.");
        } else if (memory < 1) {
            try {
                size_t total_mem = detect_available_mem(); 
                return std::floor(total_mem * memory);
            } catch (const std::exception&) {
                throw Usage("Failed to detect available memory size, please specify memory size in byte via -m option");
            }
        } else {
            return std::floor(memory);
        }
    }
};

void remove_if_exists(const std::string& file_name) {
    if (std::ifstream(file_name)) {  // Check if file exists
        if (std::remove(file_name.c_str()) == 0) {
        } else {
            std::perror("Error deleting the file");
        }
    }
}

void setup_output_dir(Options& opt) {
    if (opt.out_dir.empty()) {
        opt.out_dir = abspath("./megahit_out");
    }
    
    if (opt.mpienv_.rank == 0) {
        if (!opt.force_overwrite && !opt.test_mode && fs::exists(opt.out_dir)) {
            std::cerr << "Output directory " + opt.out_dir +
            " already exists, please change the parameter -o to another value to avoid overwriting." << std::endl;
            exit(1);
        }
    }

    if (opt.temp_dir.empty()) {
        opt.temp_dir = opt.out_dir + "/tmp";
    } else {
        try {
            // Create a temporary directory using a prefix (similar to tempfile.mkdtemp in Python)
            opt.temp_dir = fs::temp_directory_path() / ("megahit_tmp_" + std::to_string(rand()));
        } catch (const std::exception& e) {
            std::cerr << "Error creating temp directory: " << e.what() << std::endl;
            throw Usage("Failed to create temporary directory.");
        }
    }

    mkdir_if_not_exists(opt.out_dir);
    mkdir_if_not_exists(opt.temp_dir);
    mkdir_if_not_exists(opt.contig_dir());
}

void create_library_file(Options &opt) {
    FILE* lib = fopen(opt.read_lib_path().c_str(), "w");
    if (!lib) {
        perror("Failed to open read library file");
        return;
    }

    // Write PE12 entries
    for (size_t i = 0; i < opt.pe12.size(); ++i) {
        fprintf(lib, "%s\n", opt.pe12[i].c_str());
        std::string command = inpipe_cmd(opt.pe12[i]);
        if (!command.empty()) {
            fprintf(lib, "interleaved %s/inpipe.pe12.%zu\n", opt.temp_dir.c_str(), i);
        } else {
            fprintf(lib, "interleaved %s\n", opt.pe12[i].c_str());
        }
    }

    // Write PE1/PE2 entries
    for (size_t i = 0; i < opt.pe1.size(); ++i) {
        fprintf(lib, "%s,%s\n", opt.pe1[i].c_str(), opt.pe2[i].c_str());
        std::string command1 = inpipe_cmd(opt.pe1[i]);
        std::string command2 = inpipe_cmd(opt.pe2[i]);
        if (!command1.empty()) {
            fprintf(lib, "pe %s/inpipe.pe1.%zu %s/inpipe.pe2.%zu\n", opt.temp_dir.c_str(), i, opt.temp_dir.c_str(), i);
        } else {
            fprintf(lib, "pe %s %s\n", opt.pe1[i].c_str(), opt.pe2[i].c_str());
        }
    }

    // Write SE entries
    for (size_t i = 0; i < opt.se.size(); ++i) {
        fprintf(lib, "%s\n", opt.se[i].c_str());
        std::string command = inpipe_cmd(opt.se[i]);
        if (!command.empty()) {
            fprintf(lib, "se %s/inpipe.se.%zu\n", opt.temp_dir.c_str(), i);
        } else {
            fprintf(lib, "se %s\n", opt.se[i].c_str());
        }
    }

    fclose(lib);
}

void build_library(Options &opt) {
    std::vector<std::string> fifos;
    std::vector<pid_t> pipes;

    std::vector<MPI_Comm> intercomms;
    auto create_fifo = [&](const std::string& read_type, int num, int command) {
        std::string fifo_path = opt.temp_dir + "/" + "inpipe." + read_type + "." + std::to_string(num);
        remove_if_exists(fifo_path);
        //mkfifo(fifo_path.c_str(), 0644);
        fifos.push_back(fifo_path);
        
        // 创建子进程的命令和参数
        //const char* cmd;
        //if (command == 3) {
            //    cmd = "gzip";
            //} else {
                //    cmd = "bzip2";
                //}
                
                std::string cmd = abspath("./main_child"); // 子进程可执行文件的路径
                std::string wdir = abspath("./"); // 子进程可执行文件的路径
                std::cout << abspath("./main_child") << std::endl;
                //char* argv[] = {const_cast<char*>(cmd), const_cast<char*>("-cd"), const_cast<char*>(opt.se[num].c_str()), const_cast<char*>(">"), const_cast<char*>(fifo_path.c_str()), nullptr};
        char* argv[] = {const_cast<char*>(cmd.c_str()), const_cast<char*>(opt.se[num].c_str()), const_cast<char*>(fifo_path.c_str()), nullptr};
        // 使用 MPI_Comm_spawn 创建子进程
        MPI_Comm intercomm;
        MPI_Info info;
        MPI_Info_create(&info);
        MPI_Info_set(info, "wdir", wdir.c_str());
        int err = MPI_Comm_spawn("main_child", argv, 1, info, 0, MPI_COMM_SELF, &intercomm, MPI_ERRCODES_IGNORE);
        if (err != MPI_SUCCESS) {
            fprintf(stderr, "MPI_Comm_spawn failed\n");
            exit(1);
        }
        intercomms.push_back(intercomm);
    };
    /*
    auto create_fifo = [&](const std::string& read_type, int num, int &command) {
        std::string fifo_path = opt.temp_dir + "/" + "inpipe." + read_type + "." + std::to_string(num);
        remove_if_exists(fifo_path);
        fifos.push_back(fifo_path);
        
        pid_t pid = fork();
        if (pid == 0) {
            // 打开输出文件
            int fd = open(fifo_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
            if (fd == -1) {
                perror("open failed");
                exit(1);
            }
            
            // 将标准输出重定向到文件
            if (dup2(fd, STDOUT_FILENO) == -1) {
                perror("dup2 failed");
                exit(1);
            }
            close(fd);
            
            // Child process
            if (command = 3) {
                execlp("gzip", "gzip", "-cd", opt.se[num].c_str(), (char*)nullptr);
            } else {
                execlp("bzip2", "bzip2", "-cd", opt.se[num].c_str(), (char*)nullptr);
            }
            
            perror("execlp failed");
            exit(1);
        } else {
            // Parent process
            pipes.push_back(pid);
        }
    };
    */

    try {
        // Create inpipe for pe12
        for (int i = 0; i < opt.pe12.size(); ++i) {
            int cmd = execlp_cmd(opt.pe12[i]);
            if (cmd) {
                create_fifo("pe12", i, cmd);
            }
        }

        // Create inpipe for pe1
        for (int i = 0; i < opt.pe1.size(); ++i) {
            int cmd = execlp_cmd(opt.pe1[i]);
            if (cmd) {
                create_fifo("pe1", i, cmd);
            }
        }

        // Create inpipe for pe2
        for (int i = 0; i < opt.pe2.size(); ++i) {
            int cmd = execlp_cmd(opt.pe2[i]);
            if (cmd) {
                create_fifo("pe2", i, cmd);
            }
        }

        // Create inpipe for se
        for (int i = 0; i < opt.se.size(); ++i) {
            int cmd = execlp_cmd(opt.se[i]);
            if (cmd) {
                create_fifo("se", i, cmd);
            }
        }
        
        /*
        // Wait for all child processes to finish
        for (pid_t pid : pipes) {
            int status;
            waitpid(pid, &status, 0);
            if (status != 0) {
                std::cerr << "Error occurred when reading inputs\n";
                exit(status);
            }
        }
        */

        for (MPI_Comm intercomm : intercomms) {
            // 在父进程和子进程之间执行 Barrier，同步进程
            MPI_Barrier(intercomm);
            // 断开父进程与子进程之间的通信
            MPI_Comm_disconnect(&intercomm);
        }
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << '\n';
        exit(1);
    }

    std::vector<std::string> args = {"buildlib", opt.read_lib_path(), opt.read_lib_path()};

    // 直接使用 vector 来构造 argv
    std::vector<const char*> argv_vec;
    for (const auto& arg : args) {
        argv_vec.push_back(arg.c_str());  // 使用 c_str() 获取 C 风格字符串
    }

    // 传递给 main_build_lib
    main_build_lib(argv_vec.size(), const_cast<char**>(argv_vec.data()), opt.mpienv_);

    // Clean up
    for (const std::string& fifo : fifos) {
        remove_if_exists(fifo);
    }
}

void parse_option(int argc, char* argv[], Options &opt) {
    if (argc == 1) {
        std::cerr << usage_message << std::endl;
        exit(1);
    }

    for (int i = 1; i < argc; ++i) {
        std::string option = argv[i];
        if (option == "-h" || option == "--help") {
            std::cout << usage_message << std::endl;
            exit(0);
        }
        else if (option == "-o" || option == "--out-dir") {
            if (i + 1 < argc) {
                opt.out_dir = fs::absolute(argv[++i]).lexically_normal().string();
            }
        }
        else if (option == "-m" || option == "--memory") {
            if (i + 1 < argc) {
                opt.memory = std::stof(argv[++i]);
            }
        }
        else if (option == "--min-contig-len") {
            if (i + 1 < argc) {
                opt.min_contig_len = std::stoi(argv[++i]);
            }
        }
        else if (option == "-t" || option == "--num-cpu-threads") {
            if (i + 1 < argc) {
                opt.num_cpu_threads = std::stoi(argv[++i]);
            }
        }
        else if (option == "-n" || option == "--num-processes") {
            if (i + 1 < argc) {
                opt.num_processes = std::stoi(argv[++i]);
            }
        }
        else if (option == "--kmin-1pass") {
            opt.kmin_1pass = true;
        }
        else if (option == "--k-min") {
            if (i + 1 < argc) {
                opt.k_min = std::stoi(argv[++i]);
                opt.set_list_by_min_max_step = true;
                opt.auto_k = false;
            }
        }
        else if (option == "--k-max") {
            if (i + 1 < argc) {
                opt.k_max = std::stoi(argv[++i]);
                opt.set_list_by_min_max_step = true;
                opt.auto_k = false;
            }
        }
        else if (option == "--k-step") {
            if (i + 1 < argc) {
                opt.k_step = std::stoi(argv[++i]);
                opt.set_list_by_min_max_step = true;
                opt.auto_k = false;
            }
        }
        else if (option == "--k-list") {
            if (i + 1 < argc) {
                std::string value = argv[++i];
                std::stringstream ss(value);
                std::string temp;
                while (std::getline(ss, temp, ',')) {
                    opt.k_list.push_back(std::stoi(temp));
                }
                std::sort(opt.k_list.begin(), opt.k_list.end());
                opt.auto_k = false;
                opt.set_list_by_min_max_step = false;
            }
        }
        else if (option == "--min-count") {
            if (i + 1 < argc) {
                opt.min_count = std::stoi(argv[++i]);
            }
        }
        else if (option == "--max-tip-len") {
            if (i + 1 < argc) {
                opt.max_tip_len = std::stoi(argv[++i]);
            }
        }
        else if (option == "--merge-level") {
            if (i + 1 < argc) {
                std::string value = argv[++i];
                std::stringstream ss(value);
                std::string temp;
                std::getline(ss, temp, ',');
                opt.merge_len = std::stoi(temp);
                std::getline(ss, temp, ',');
                opt.merge_similar = std::stof(temp);
            }
        }
        else if (option == "--prune-level") {
            if (i + 1 < argc) {
                opt.prune_level = std::stoi(argv[++i]);
            }
        }
        else if (option == "--prune-depth") {
            if (i + 1 < argc) {
                opt.prune_depth = std::stof(argv[++i]);
            }
        }
        else if (option == "--bubble-level") {
            if (i + 1 < argc) {
                opt.bubble_level = std::stoi(argv[++i]);
            }
        }
        else if (option == "--no-mercy") {
            opt.no_mercy = true;
        }
        else if (option == "--no-local") {
            opt.no_local = true;
        }
        else if (option == "--disconnect-ratio") {
            if (i + 1 < argc) {
                opt.disconnect_ratio = std::stof(argv[++i]);
            }
        }
        else if (option == "--low-local-ratio") {
            if (i + 1 < argc) {
                opt.low_local_ratio = std::stof(argv[++i]);
            }
        }
        else if (option == "--cleaning-rounds") {
            if (i + 1 < argc) {
                opt.cleaning_rounds = std::stoi(argv[++i]);
            }
        }
        else if (option == "--keep-tmp-files") {
            opt.keep_tmp_files = true;
        }
        else if (option == "--mem-flag") {
            if (i + 1 < argc) {
                opt.mem_flag = std::stoi(argv[++i]);
            }
        }
        else if (option == "-v" || option == "--version") {
            std::cout << "Version info goes here" << std::endl;
            exit(0);
        }
        else if (option == "--verbose") {
            opt.verbose = true;
        }
        else if (option == "--continue") {
            opt.continue_mode = true;
        }
        else if (option == "--out-prefix") {
            if (i + 1 < argc) {
                opt.out_prefix = argv[++i];
            }
        }
        else if (option == "--tmp-dir") {
            if (i + 1 < argc) {
                opt.temp_dir = fs::absolute(argv[++i]).lexically_normal().string();
            }
        }
        else if (option == "--force") {
            opt.force_overwrite = true;
        }
        else if (option == "--test") {
            opt.test_mode = true;
        }
        else if (option == "--no-hw-accel") {
            opt.hw_accel = false;
            opt.has_popcnt = false;
        }
        else if (option == "-r" || option == "--read") {
            if (i + 1 < argc) {
                std::string value = argv[++i];
                std::stringstream ss(value);
                std::string temp;
                while (std::getline(ss, temp, ',')) {
                    opt.se.push_back(fs::absolute(temp).lexically_normal().string());
                }
            }
        }
        else if (option == "-1") {
            if (i + 1 < argc) {
                std::string value = argv[++i];
                std::stringstream ss(value);
                std::string temp;
                while (std::getline(ss, temp, ',')) {
                    opt.pe1.push_back(fs::absolute(temp).lexically_normal().string());
                }
            }
        }
        else if (option == "-2") {
            if (i + 1 < argc) {
                std::string value = argv[++i];
                std::stringstream ss(value);
                std::string temp;
                while (std::getline(ss, temp, ',')) {
                    opt.pe2.push_back(fs::absolute(temp).lexically_normal().string());
                }
            }
        }
        else if (option == "--12") {
            if (i + 1 < argc) {
                std::string value = argv[++i];
                std::stringstream ss(value);
                std::string temp;
                while (std::getline(ss, temp, ',')) {
                    opt.pe12.push_back(fs::absolute(temp).lexically_normal().string());
                }
            }
        }
        else if (option == "--presets") {
            if (i + 1 < argc) {
                opt.presets = argv[++i];
            }
        }
        else {
            std::cerr << "Invalid option " << option << std::endl;
            exit(1);
        }
    }
}

void check_and_correct_option(Options& opt) {
    // set mode
    if (opt.memory < 0) {
        throw Usage("-m cannot be less than 0!");
    }

    if (opt.presets != "") {
        opt.auto_k = true;

        if (opt.presets == "meta-sensitive") {
            opt.min_count = 1;
            opt.k_list = {21, 29, 39, 49, 59, 69, 79, 89, 99, 109, 119, 129, 141};
            opt.set_list_by_min_max_step = false;
        }
        else if (opt.presets == "meta-large") {
            opt.min_count = 1;
            opt.k_min = 27;
            opt.k_max = 127;
            opt.k_step = 10;
            opt.set_list_by_min_max_step = true;
        }
        else {
            throw Usage("Invalid preset: " + opt.presets);
        }
    }

    if (opt.set_list_by_min_max_step) {
        if (opt.k_step % 2 == 1) {
            throw Usage("k-step must be even number!");
        }
        if (opt.k_min > opt.k_max) {
            throw Usage("Error: k_min > k_max!");
        }

        opt.k_list.clear();
        int k = opt.k_min;
        while (k < opt.k_max) {
            opt.k_list.push_back(k);
            k += opt.k_step;
        }
        opt.k_list.push_back(opt.k_max);
    }

    if (opt.k_list.empty()) {
        throw Usage("k list should not be empty!");
    }

    if (opt.k_list.front() < 15 || opt.k_list.back() > 255) { // example max_k_allowed
        throw Usage("All k's should be in range [15, 255]");
    }

    for (int k : opt.k_list) {
        if (k % 2 == 0) {
            throw Usage("All k must be odd number!");
        }
    }

    for (size_t i = 1; i < opt.k_list.size(); ++i) {
        if (opt.k_list[i] - opt.k_list[i - 1] > 28) {
            throw Usage("--k-step (adjacent k difference) must be <= 28");
        }
    }

    opt.k_min = opt.k_list.front();
    opt.k_max = opt.k_list.back();

    if (opt.k_max < opt.k_min) {
        throw Usage("--k-min should not be larger than --k-max.");
    }
    if (opt.min_count <= 0) {
        throw Usage("--min-count must be greater than 0.");
    }
    else if (opt.min_count == 1) {
        opt.kmin_1pass = true;
        opt.no_mercy = true;
    }
    if (opt.prune_level < 0 || opt.prune_level > 3) {
        throw Usage("--prune-level must be in 0-3.");
    }
    if (opt.merge_len < 0) {
        throw Usage("--merge-level: length must be >= 0");
    }
    if (opt.merge_similar < 0 || opt.merge_similar > 1) {
        throw Usage("--merge-level: similarity must be in [0, 1]");
    }
    if (opt.disconnect_ratio < 0 || opt.disconnect_ratio > 0.5) {
        throw Usage("--disconnect-ratio should be in [0, 0.5].");
    }
    if (opt.low_local_ratio <= 0 || opt.low_local_ratio > 0.5) {
        throw Usage("--low-local-ratio should be in (0, 0.5].");
    }
    if (opt.cleaning_rounds <= 0) {
        throw Usage("--cleaning-rounds must be >= 1");
    }

    // Assuming we use `sysconf` to get the number of CPU cores instead of `os.sched_getaffinity`
    long cpu_cores = sysconf(_SC_NPROCESSORS_ONLN);
    if (opt.num_cpu_threads > cpu_cores) {
        std::cerr << "Warning: Maximum number of available CPU threads is " << cpu_cores << ".\n";
        std::cerr << "Number of threads is reset to " << cpu_cores << ".\n";
        opt.num_cpu_threads = cpu_cores;
    }
    if (opt.num_cpu_threads == 0) {
        opt.num_cpu_threads = cpu_cores;
    }

    if (opt.prune_depth < 0 && opt.prune_level < 3) {
        opt.prune_depth = opt.min_count;
    }
    if (opt.bubble_level < 0) {
        std::cerr << "Warning: Reset bubble level to 0.\n";
        opt.bubble_level = 0;
    }
    if (opt.bubble_level > 2) {
        std::cerr << "Warning: Reset bubble level to 2.\n";
        opt.bubble_level = 2;
    }
}

void build_first_graph(Options& opt) {
    if (!opt.kmin_1pass) {

        std::vector<std::string> args_count = {"count", "-k", std::to_string(opt.k_min)
                                                , "-m", std::to_string(opt.min_count)
                                                , "--host_mem", std::to_string(opt.host_mem())
                                                , "--mem_flag", std::to_string(opt.mem_flag)
                                                , "--output_prefix", graph_prefix(opt.temp_dir, opt.k_min)
                                                , "--num_cpu_threads", std::to_string(opt.num_cpu_threads)
                                                , "--read_lib_file", opt.read_lib_path()};
        
        std::vector<const char*> argv_vec;
        for (const auto& arg : args_count) {
            argv_vec.push_back(arg.c_str());
        }

        main_kmer_count(argv_vec.size(), const_cast<char**>(argv_vec.data()), opt.mpienv_);

        
    } else {
        std::vector<std::string> args_1pass = {"read2sdbg"  , "-k", std::to_string(opt.k_min)
                                                            , "-m", std::to_string(opt.min_count)
                                                            , "--host_mem", std::to_string(opt.host_mem())
                                                            , "--mem_flag", std::to_string(opt.mem_flag)
                                                            , "--output_prefix", graph_prefix(opt.temp_dir, opt.k_min)
                                                            , "--num_cpu_threads", std::to_string(opt.num_cpu_threads)
                                                            , "--read_lib_file", opt.read_lib_path()};
        //no mercy TODO
        std::vector<const char*> argv_vec;
        for (const auto& arg : args_1pass) {
            argv_vec.push_back(arg.c_str());
        }

        main_read2sdbg(argv_vec.size(), const_cast<char**>(argv_vec.data()), opt.mpienv_);
    }
    


    if (!opt.keep_tmp_files) {
        //remove_temp_after_build(opt.k_min);
    }
    
}

void setup_logger(Options &opt) {
}

void assemble(Options& opt, int cur_k) {
    int min_standalone = std::max(std::min(opt.k_max * 3 - 1, (int)(opt.min_contig_len * 1.5)), opt.min_contig_len);
    if (opt.max_tip_len >= 0)
        min_standalone = std::max(opt.max_tip_len + opt.k_max - 1, opt.min_contig_len);

    std::vector<std::string> args_assemble = {"assemble"    , "-s", graph_prefix(opt.temp_dir, cur_k)
                                                            , "-o", contig_prefix(opt.contig_dir(), cur_k)
                                                            , "-t", std::to_string(opt.num_cpu_threads)
                                                            , "--min_standalone", std::to_string(min_standalone)
                                                            , "--prune_level", std::to_string(opt.prune_level)
                                                            , "--merge_len", std::to_string(opt.merge_len)
                                                            , "--merge_similar", std::to_string(opt.merge_similar)
                                                            , "--cleaning_rounds", std::to_string(opt.cleaning_rounds)
                                                            , "--disconnect_ratio", std::to_string(opt.disconnect_ratio)
                                                            , "--low_local_ratio", std::to_string(opt.low_local_ratio)
                                                            , "--cleaning_rounds", std::to_string(opt.cleaning_rounds)
                                                            , "--min_depth", std::to_string(opt.prune_depth)
                                                            , "--bubble_level", std::to_string(opt.bubble_level)};
    
    if ((opt.max_tip_len == -1) && ((cur_k * 3 - 1) > (opt.min_contig_len * 1.5))) {
        args_assemble.push_back("--max_tip_len");
        args_assemble.push_back(std::to_string(std::max(1, (int)(opt.min_contig_len * 1.5 + 1 - cur_k))));
    } else {
        args_assemble.push_back("--max_tip_len");
        args_assemble.push_back(std::to_string(opt.max_tip_len));
    }

    if (cur_k < opt.k_max)
        args_assemble.push_back("--careful_bubble");
    
    if (cur_k == opt.k_max)
        args_assemble.push_back("--is_final_round");
        
    if (opt.no_local)
        args_assemble.push_back("--output_standalone");

    
    std::vector<const char*> asm_args;
    for (const auto& arg : args_assemble) {
        asm_args.push_back(arg.c_str());
    }

    main_assemble(asm_args.size(), const_cast<char**>(asm_args.data()), opt.mpienv_);

    if (!opt.keep_tmp_files && cur_k != opt.k_max) {
        //remove_temp_after_assemble();
    }
}

int main(int argc, char **argv) {
    Options opt;
    opt.mpienv_.init(argc, argv);

    
    parse_option(argc, argv, opt);
    setup_output_dir(opt);
    setup_logger(opt);
    
    check_and_correct_option(opt);

    if (opt.mpienv_.rank == 0) {
        create_library_file(opt);
        build_library(opt);
    }

    MPI_Barrier(MPI_COMM_WORLD); // Barrier
    build_first_graph(opt);

    MPI_Barrier(MPI_COMM_WORLD); // Barrier
    assemble(opt, opt.k_min);

    opt.mpienv_.finalize();
    return 0;
}
