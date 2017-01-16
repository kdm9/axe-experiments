#!/usr/bin/awk -f
{
    if ($0 ~ /^#/) {
        print $0;
    } else {
        printf("%s\t%s\n", $2, $1);
    }
}
