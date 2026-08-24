import sys
import time

import highspy

for path in sys.argv[1:]:
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    h.readModel(path)
    t0 = time.perf_counter()
    h.run()
    dt = (time.perf_counter() - t0) * 1000.0
    info = h.getInfo()
    name = path.split("/")[-1].rsplit(".", 1)[0]
    print(f"{name},{dt:.1f},{info.objective_function_value:.10g},"
          f"{h.modelStatusToString(h.getModelStatus())}")
