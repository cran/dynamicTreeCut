if (exists("printFlush")) { remove(printFlush); collect_garbage(); }
printFlush = function(...)
{
#  x = print(...)
  cat(...); cat("\n");
  if (exists("flush.console")) flush.console();
}

if (exists("indentSpaces")) { remove(indentSpaces); collect_garbage(); }
indentSpaces = function(indent = 0)
{
  if (indent>0) 
    {
      spaces = paste(rep("  ", times=indent), collapse="");
    } else
    {
      spaces = "";
    }
  spaces;
}

