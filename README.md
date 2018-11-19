# FPWS1819
Fortgeschrittenen Praktikum WS 18/19

Hallo Alexander,

ich schreibe dir einfach mal ein paar Tips in die README.

Falls dir die Windows Umgebung ein Terminal bietet:
Die vier Hauptbefehle sind (Linux):

git pull        -> Lade die aktuelle Version herunter
git add         -> Markiere files fuer den commit
git commit      -> Erstellt einen commit (Aktuellen Stand des Repositories speichern)
git push        -> Laed alle commits hoch

Tipps:

git commit commitet nur das was mit git add geadded wurde.
git commit -m "message" adds the commit message "message" to the commit

git add -A added alles im Ordner (nicht empfohlen).
git add -u added alles, was bereits einmal commitet, aber veraendert wurde.
Den -u kannst du verwenden um veraenderte Dateien zu adden.

NIEMALS einfach ein file loeschen, bewegen, oder umbennenen, das bereits Teil des git ist.
Dafuer benutzt man git rm und git mv. (rm=remove mv=move)
In Windows mag das automatisch funktionieren.
