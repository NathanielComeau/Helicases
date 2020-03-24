# Notes on Samcomp

Ibrahim has repeatedly mentionned that we should try and get Samcomp up and running (among other things).

I downloaded the source code from here: https://sourceforge.net/projects/samcomp/

I found the samcomp paper here: https://doi.org/10.1371/journal.pone.0059190

In the Samcomp source code the README mentions how it compresses quality scores:

```text
Qualities get much the same treatment here as in fqzcomp - the
previous byte, the maximum of the previous two[1], and an overall
deviation score are used as context to predict the next value. Unlike
fqzcomp there is currently no per-position context. [1]The reason for
this is that low quality values can often occur early on just due to
abherrant calls. Eg a quality string maybe 38 38 37 37 12 37 37.
The 37s after the 12 have much more in common with the previous 37 and
not the previous 12, so using the maximum of previous 2 improves
predictions and decreases encoding costs.
```

Ibrahim potentially implemented this in Deeze: https://github.com/inumanag/deez/blob/master/Streams/SAMCompStream.h



Our goal: get it running for a baseline of quality score compression performance.



