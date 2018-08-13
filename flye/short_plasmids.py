import multiprocessing
import signal

from flye.alignment import SynchronizedSamReader
import flye.fasta_parser as fp


def _thread_worker(aln_reader, results_queue, error_queue):
    try:
        aln_reader.init_reading()

        while not aln_reader.is_eof():
            ctg_id, ctg_aln = aln_reader.get_chunk()
            if ctg_id is None:
                break

    except Exception as e:
        error_queue.put(e)


def find_unmapped_reads(alignment_path, contigs_path, num_proc, min_aln_length):
    aln_reader = SynchronizedSamReader(alignment_path,
                                       fp.read_fasta_dict(contigs_path),
                                       min_aln_length)
    manager = multiprocessing.Manager()
    results_queue = manager.Queue()
    error_queue = manager.Queue()

    orig_sigint = signal.signal(signal.SIGINT, signal.SIG_IGN)
    threads = []
    for _ in xrange(num_proc):
        threads.append(multiprocessing.Process(target=_thread_worker,
                                               args=(aln_reader,
                                                     results_queue,
                                                     error_queue)))
    signal.signal(signal.SIGINT, orig_sigint)

    for t in threads:
        t.start()
    try:
        for t in threads:
            t.join()
    except KeyboardInterrupt:
        for t in threads:
            t.terminate()

    if not error_queue.empty():
        raise error_queue.get()
