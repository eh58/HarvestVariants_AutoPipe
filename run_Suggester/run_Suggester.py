import os
import sys
import tempfile
import subprocess

def is_warning_line(line):
    return line.startswith("WARNING")

# Modify the argument check and help message
if len(sys.argv) < 3:
    print("Usage: {} <path_to_jar> <input_xml_file> [output_directory]".format(sys.argv[0]))
    sys.exit(1)

path_to_jar = sys.argv[1]
input_xml_file = sys.argv[2]

# # Set the output_directory to the provided value or use the default
# if len(sys.argv) == 4:
#     output_directory = sys.argv[3]
# else:
#     output_directory = "./Data/suggester/"

output_directory = "./Data/suggester/"
os.makedirs(output_directory, exist_ok=True)

log4j2_config = tempfile.NamedTemporaryFile(delete=False)
log4j2_config.write(b'''<?xml version="1.0" encoding="UTF-8"?>
<Configuration status="WARN">
  <Appenders>
    <Console name="Console" target="SYSTEM_ERR">
      <PatternLayout pattern="%d{HH:mm:ss.SSS} [%t] %-5level %logger{36} - %msg%n"/>
    </Console>
  </Appenders>
  <Loggers>
    <Root level="error">
      <AppenderRef ref="Console"/>
    </Root>
  </Loggers>
</Configuration>
''')
log4j2_config.close()

if input_xml_file.endswith("_pp.xml"):
    file_identifier = os.path.basename(input_xml_file)[:-7]

    output_file = os.path.join(output_directory, "{}.txt".format(file_identifier))
    log_file = os.path.join(output_directory, "{}.log".format(file_identifier))

    cmd = [
        "java",
        "-Dlog4j.configurationFile={}".format(log4j2_config.name),
        "-jar",
        path_to_jar,
        input_xml_file,
    ]

    # Run the command and capture the standard output and standard error
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Write non-warning lines from standard output to the output file
    with open(output_file, "w") as stdout_file:
        for line in result.stdout.splitlines():
            if not is_warning_line(line):
                stdout_file.write(line + '\n')

    # Write non-warning lines from standard error to the log file
    with open(log_file, "w") as stderr_file:
        for line in result.stderr.splitlines():
            if not is_warning_line(line):
                stderr_file.write(line + '\n')

os.remove(log4j2_config.name)
