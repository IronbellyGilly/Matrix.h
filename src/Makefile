CC = gcc
CFLAGS = -Wall -Werror -Wextra
TEST_LIBS = -lcheck -lm -lpthread
OBJECTS=s21_matrix.o
SOURSES=s21_matrix.c
GCOVFLAGS = -fprofile-arcs -ftest-coverage

all: clean s21_matrix.a 

s21_matrix.a: s21_matrix.o
	$(CC) -c $(SOURSES)
	ar rcs s21_matrix.a $(OBJECTS)

test: test.c s21_matrix.a
	$(CC) $(CFLAGS) -c test.c -o test.o
	$(CC) test.o s21_matrix.a $(TEST_LIBS) -o test && ./test 

gcov_report: s21_matrix.a 
	$(CC) $(CFLAGS) --coverage test.c s21_matrix.c s21_matrix.a $(TEST_LIBS) -o gcov_test
	./gcov_test
	lcov -t "gcov_test" -o gcov_test.info -c -d .
	genhtml -o report/ gcov_test.info
	open ./report/index.html

test_val:
	$(CC) *.c -o run $(TEST_LIBS) -lgcov --coverage
	sh valgrind.sh ./run

clean:
	rm -rf *.o *.a *.a.out gcov_test ./GcovReport *.gcno *gcda report/ *.info *.dSYM ./CPPLINT* test

rebuild: clean all

leaks: test
	CK_FORK=no leaks --atExit -- ./test

clang: test s21_matrix.a
	cp ../materials/linters/.clang-format .clang-format
	clang-format -n *.c *.h
	clang-format -n test.c
	clang-format -i *.c *.h
	clang-format -i test*.c