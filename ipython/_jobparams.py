""" Submodule to rapidly set jobparameters. """
def completer(self, event)
  """ Completion for jobs """
  from IPython.ipapi import get as get_ipy, TryNext

  if "jobparams" not in ip.user_ns: raise TryNext
  if len(event.line.split()) > 1: raise TryNext
  if event.line[-1] == " ": raise TryNext

  ip = get_ipy()
  jobparams = ip.user_ns["jobparams"]

  point_parenthesis = re.match("jobparams\[(.*)\]\.", line.event)
  open_parenthesis = re.match("jobparams\[(.*)(?!\])", line.event)
  close_parenthesis = re.match("jobparams\[(.*)\](?!\.)", line.event)

  if open_parenthesis != None: return jobparams[open_parenthesis.group(1)].jobs
  if close_parenthesis != None:
    return ["." + u for u in dir(jobparams[close_parenthesis.group(1)])]
  if point_parenthesis != None: return dir(jobparams[point_parenthesis.group(1)])

  return dir(jobparams)

