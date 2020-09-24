import threading
import time
from queue import Queue

print_lock = threading.Lock()


class MyThread(threading.Thread):
    def __init__(self, input_queue, output_queue):
        super().__init__()
        self.queue = input_queue
        self.output_queue = output_queue
        self.daemon = True

    def run(self):
        while True:
            val = self.queue.get()
            if val is None:  # If you send `None`, the thread will exit.
                return
            self.do_thing_with_message(val)

    def do_thing_with_message(self, message):
        with print_lock:
            print(threading.currentThread().getName(), "Received {}".format(message))
            self.output_queue.put(message)


if __name__ == '__main__':
    threads = []
    for t in range(10):
        q = Queue()
        output_queue = Queue()
        threads.append(MyThread(q, output_queue))
        threads[t].start()
        time.sleep(0.1)

    for x in range(1000):
        for t in threads:
            t.queue.put(x)

    for t in threads:
        while not t.output_queue.empty():
            data = t.output_queue.get()
            print(data)
        t.queue.put(None)
        t.join()
