
#include <iostream>
#include <cassert>
#include <process.h>
#include <Windows.h>

class Runnable {
public:
	virtual void run() = 0;
	virtual ~Runnable() = 0;
};

Runnable::~Runnable(){};

class Thread :public Runnable {
public:
	Thread();
	Thread(Runnable*);
	virtual ~Thread();
	bool isAlive();
	void join();
	void start();
	Thread& operator=(const Thread&);
private:
	HANDLE hThread;
	unsigned wThreadID;
	Runnable *m_runnable;
	bool *m_isAlive;
	void* result;
	virtual void run(){}
	void setCompleted();
	void setIsAlive(bool);
	static unsigned WINAPI startThreadRunnable(LPVOID pVoid);
	static unsigned WINAPI startThread(LPVOID pVoid);
};

Thread::Thread() : m_runnable(NULL) {
	m_isAlive = new bool();
	setIsAlive(false);
	hThread = (HANDLE)_beginthreadex(NULL,0,Thread::startThread, (LPVOID)this, CREATE_SUSPENDED, &wThreadID);
}

Thread::Thread(Runnable* runnable) :m_runnable(runnable) {
	m_isAlive = new bool();
	setIsAlive(false);
	hThread = (HANDLE)_beginthreadex(NULL,0,Thread::startThreadRunnable, (LPVOID)this, CREATE_SUSPENDED, &wThreadID);
}

Thread::~Thread() {
	if(wThreadID != GetCurrentThreadId()) {
		CloseHandle(hThread);
	}
}

unsigned WINAPI Thread::startThread(LPVOID pVoid) {
	Thread* aThread = static_cast<Thread*>(pVoid);
	aThread->run();
	aThread->setCompleted();
	return 0;
}

unsigned WINAPI Thread::startThreadRunnable(LPVOID pVoid) {
	Thread* runnableThread = static_cast<Thread*>(pVoid);
	runnableThread->m_runnable->run();
	runnableThread->setCompleted();
	return reinterpret_cast<unsigned>(runnableThread->result);
}
bool Thread::isAlive() {
	return *m_isAlive;
}

void Thread::join() {
	while(isAlive())
		WaitForMultipleObjects(1, &hThread, TRUE, 0);
}

void Thread::setCompleted() {
	setIsAlive(false);
}

void Thread::setIsAlive(bool isAlive) {
	*m_isAlive = isAlive;
}
void Thread::start() {
	setIsAlive(true);
	assert(hThread);
	ResumeThread(hThread);
}

Thread& Thread::operator=(const Thread& t) {
	if(this == &t) {
		return *this;
	}
	hThread = t.hThread;
	m_runnable = t.m_runnable;
	m_isAlive = t.m_isAlive;
	return *this;
}
