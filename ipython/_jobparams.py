""" Submodule to rapidly set jobparameters. """
def completer(self, event):
  """ Completion for jobs """
  from IPython.ipapi import get as get_ipy, TryNext
  from os.path import join
  from re import match

  ip = get_ipy()
  if "jobparams" not in ip.user_ns: raise TryNext
  if len(event.line.split()) > 1: raise TryNext
  if event.line[-1] == " ": raise TryNext

  jobparams = ip.user_ns["jobparams"]

  open_parenthesis = match("jobparams\[(\"|\')(.*)(?!\"|\')", event.line)
  close_parenthesis = match("jobparams\[(\"|\')(.*)(?:\"|\')", event.line)

  if open_parenthesis != None:
    name = open_parenthesis.group(2)
    if len(name) == 0:  return [event.line + u for u in jobparams.jobdict.children.keys()]
    first = ""
    symbol = open_parenthesis.group(1)
    if name[0] == '/': first, name = name[0], name[1:]
    if '/' in name:
      if name[-1] != '/': name = name[:-name[::-1].find('/')]
      if (first + name) not in jobparams.jobdict: return []
      jobdict = jobparams.jobdict[first+name]
      if len(jobdict.children) == 0:
        return [event.line.replace(open_parenthesis.group(2), first+name+symbol+"]")]
      return [ event.line.replace(open_parenthesis.group(2), join(first+name, u)) + "/" \
               for u in jobdict.children.keys() ]
    return ["jobparams[\"" + u for u in jobparams.jobdict.children.keys()]
  if close_parenthesis != None:
    symbol = close_parenthesis.group(1)
    start = "jobparams[" + symbol + close_parenthesis.group(2) + symbol + "]."
    j = jobparams[close_parenthesis.group(1)]
    return [start + u for u in dir(j)]

  if event.line[-1] == '[':
    return ["jobparams[\'" + u + "/" for u in jobparams.jobdict.children.keys()]
  result = [ "jobparams." + u for u in dir(jobparams) ]
  if len(jobparams.jobdict.children) > 0: result.append("jobparams[\"")
  return result

