# Linear-algebra
## Dependencies: oneAPI

### 下载安装

[oneAPI官网](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html#gs.852ahl)

进入官网后注册账户即可进行下载，这里我选择的是命令行下载、本地安装

```shell
wget https://registrationcenter-download.intel.com/akdlm/irc_nas/17977/l_BaseKit_p_2021.3.0.3219_offline.sh

sudo bash l_BaseKit_p_2021.3.0.3219_offline.sh
# 这里我选择安装在根目录下所以加上了sudo，安装在home目录下应该也是可以的，不过要记住路径以便后续环境变量配置
# 安装过程全都默认即可，也可以自己更换安装目录
```

### 关于目录

![安装目录](https://raw.githubusercontent.com/Chen-WH/PicGo/main/Typora/202109301536568.png)

- oneapi
  - mkl/
    - 2021.3.0/
      - include/
      - lib/
      - 其他
    - latest/ --> 2021.3.0/
  - mpi/
  - tbb/
  - 其他工具包以及安装程序等

以mkl工具包为例，底下存放有mkl包的各个历史版本，在我的截图这里可以看到只有 2021.3.0 一个版本和 latest 链接，latest链接随更新会保持链接到最新版本，因此在工程文件链接库目录的时候可以选择 latest 路径，当然指定某一版本也是可以的。

### 环境配置

```shell
# setvars.sh 脚本通过在各自的 oneAPI 文件夹中找到每个 <install-dir>/latest/env/vars.sh 脚本来设置用于 oneAPI 工具包的环境变量。
source /opt/intel/oneapi/setvars.sh intel64
# 为了避免每次打开一个新的终端都需要重新运行脚本，可以在 ${HOME}/.bashrc中添加一句 source /opt/intel/oneapi/setvars.sh intel64
# 运行这句后oneAPI会添加所有工具包，可能会污染你的工作环境，比如我就遇到了 intelpython 影响 rospy 的问题，可以写一个 config.txt 文件自定义加载的环境
sudo vi /opt/intel/oneapi/config.txt # 在指定目录下新建了一个配置文件，内容如下
intelpython=exclude # 这句排除了intelpython，可以自己进行设置，这里只是举了我的例子
# 那么source语句需要进行适当修改如下
source /opt/intel/oneapi/setvars.sh --config="/opt/intel/oneapi/config.txt"

## 其他例子如下
mkl = 1.1 # 指定使用的mkl版本，如果不指定默认使用最新版本
default = exclude # 指定默认为排除所有，如果不设置的话默认添加所有工具包

## 去除命令行echo输出
# source setvars.sh后每次命令行会出现很多反馈，如果像我一样强迫症不希望显示的话可以这样修改source语句
source /opt/intel/oneapi/setvars.sh --config="/opt/intel/oneapi/config.txt" > /dev/null
# /dev/null 是类Unix系统中的一个特殊文件设备，他的作用是接受一切输入它的数据并丢弃这些数据。通常被当做垃圾桶来用。
```

### 如何使用

#### 命令行编译

```c
#include "mkl.h"
// 首先include头文件，具体编译选项非常复杂，参考下方网站进行补充
```

[编译选项参考](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl/link-line-advisor.html)

#### 使用Cmake

```cmake
# CMakeLists.txt
cmake_minimum_required(VERSION 3.0.2)
project(eigen_mkl)
# set(CMAKE_BUILD_TYPE "Release" )
# set(CMAKE_BUILD_TYPE "Debug" )
set(CMAKE_CXX_FLAGS "-O3" )

# include 头文件
include_directories
  ${catkin_INCLUDE_DIRS}
  /opt/intel/oneapi/mkl/latest/include/
)

# 链接库目录
link_directories(
  ${catkin_LIB_DIRS}
  /opt/intel/oneapi/mkl/latest/lib/
 )
 
## 文件名是eigen_mkl.cpp
add_executable(Matrix_mkl src/eigen_mkl.cpp)
# 链接具体的库 libmkl_rt，注意.so是动态库，当然也可以选择静态库，不过文件会大些但更稳定。
target_link_libraries(eigen_mkl
libmkl_rt.so
)
```

## 
