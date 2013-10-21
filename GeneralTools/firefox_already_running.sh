#!/bin/sh

LOCK=$(find $HOME/.mozilla/firefox/ -name lock)
PLOCK=$(find $HOME/.mozilla/firefox/ -name \.parentlock)

echo "Deleting Lock file:	$LOCK"
rm -f $LOCK
echo "Deleting Parent Lock file: $PLOCK"
rm -f $PLOCK

echo "To see why these files had to be deleted, see: http://www.mattcutts.com/blog/how-to-fix-firefox-is-already-running-error/"
echo "Try running firefox again..."
