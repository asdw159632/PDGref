cmake_minimum_required(VERSION 3.10)

# 项目名称和版本
	project(PDGref VERSION 1.0)

# 找到 ROOT 包
	find_package(ROOT REQUIRED)

	if(ROOT_FOUND)
		message(STATUS "ROOT version ${ROOT_VERSION}")

		# 获取 ROOT 库cflags
		execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CFLAGS
			OUTPUT_STRIP_TRAILING_WHITESPACE)

		# 获取 ROOT 库libs
		execute_process(COMMAND root-config --libs OUTPUT_VARIABLE ROOT_LIBS
			OUTPUT_STRIP_TRAILING_WHITESPACE)

	else()
		message(FATAL_ERROR
			" There is no CERN ROOT found, please check your ROOT
			support\n"
			" Program must be built with ROOT")
	endif()

	# 将 ROOT 编译标志添加到编译器标志
	set(CMAKE_CXX_FLAGS "${ROOT_CFLAGS}")

	# 添加可执行文件
	add_executable(PDGref src/main.cxx)

	set(CMAKE_CXX_STANDARD_LIBRARIES "-lstdc++")

	# 链接 ROOT 库
	target_link_libraries(PDGref ${ROOT_LIBS} -lEG)

