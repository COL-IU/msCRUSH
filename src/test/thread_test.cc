#include <iostream>
#include <thread>
using namespace std;

void call_from_thread(int *n) {
  cout << "Hello world!" << endl;
  cout << *n << endl;
  (*n)++;
}
int main () {
  int num = 0;
  cout << "num: " << num << endl;
  thread t1(call_from_thread, &num);
  thread t2(call_from_thread, &num);
  thread t3(call_from_thread, &num);
  t1.join();
  t2.join();
  t3.join();
  cout << "num: " << num << endl;
  return 0;
}
