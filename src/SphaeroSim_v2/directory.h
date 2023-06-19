#ifndef __DIRECTORY_H
#define __DIRECTORY_H

#include <string>
#include <vector>

#ifndef WIN32
#ifndef __linux__
#define __linux__
#endif
#endif

#ifdef __linux__

#include <sys/types.h>
#include <dirent.h>

#else

#include <windows.h>
#include <tchar.h>

#endif // __linux__

#ifndef TCHAR
#define TCHAR char
#endif

#ifndef _T
#define _T
#endif


//! Stellt eine OS unabhängige Klasse dar für Ordnerzugriff
/*!
  Mit dieser Klasse können Elemente aus Ordnern
  ausgelesen werden. Mit der Compilereinstellung
  in defines.h kann das OS angebenen werden.

  Diese Klasse iteriert über die Elemente aus einem
  Ordner, dazu sind die Standartiteratorfunktionen
  überladen.

  \author Christian Dernehl
  \date 2008
 */
class Directory {
public:

    //! Konstruktor auf kein Verzeichnis
    Directory();
    //! Konstruktor mit Verzeichnisinitialisierung
    /*!
      Konstrukiert ein Objekt und versucht dieses
      mit dem angebenen Pfad zu öffnen.

      @param[in] path Der Pfad zum Verzeichnis

      \see open_dir(const char* path)

      \warning Der Pfad muss gültig sein sonst ist das Verhalten nicht definiert
     */
    Directory(const TCHAR *path);

    //! Copy Konstruktor
    Directory(const Directory &rhs);

    //! Destruktor
    ~Directory();

    //! Zuweisungsopertor
    /*!
      Führt eine ähnliche Funktion wie der
      Copy Konstruktor aus.
     */
    Directory &operator=(const Directory &rhs);
    //! Präfix Inkremtationsoperator

    //! öffnet ein Verzeichnis
    /*!
      Wird ein ungültiger Pfad angegeben oder
      kann der Verzeichnispointer nicht vom
      OS bereitgestellt werden wird ein wahr
      Wert zurückgegeben

      @param[in] path Der Pfad zum Verzeichnis

      \return true wenn erfolgreich, false wenn Fehler
     */
    bool open_dir(const TCHAR *path);

#ifdef __linux__
    //! Alle Verzeichniselemente
    std::vector<dirent> elements;
    //! Das aktuelle Element
    std::vector<dirent>::iterator it;
#else
    //! Alle Verzeichniselemente
    std::vector<WIN32_FIND_DATA> elements;
#endif

    // return the filename for an element
    std::string getFileNameFromIndex(int iIndex);
    //! Gibt den Pfad gekoppelt mit dem Verzeichnis zurück
    /*!
      Verknüpft Verzeichnispfad mit dem Pfad des aktuellen
      Elementes.

      \warning Diese Funktion ist noch nicht für alle Fälle abgedeckt
     */
    std::string get_file_full_path();

    // returns true if the element iIndex is a folder
    bool isFolder(int iIndex);

    // returns true if pathname exists as a folder
    static bool folderExists(const std::string &strPathname);

    //! Sortiert die Dateinamen alphanumerisch
    void sort();

    //! Sortiert die Dateinamen alphanumerisch, dabei stehen kuerzere Dateinamen immer weiter vorne
    void sort_short_filenames_first();

    //! überprüft ob eine Datei existiert
    /*!
      öffnet eine Filestream auf die Datei und testet ob die Datei
      gelesen werden kann.

      @param[in] filename Der Pfad zur Datei
      \return True wenn die Datei existiert und geöffnet werden kann
     */
    static bool file_exists(const std::string &filename);

    //! Erstellt ein neues Verzeichnis
    /*!
      Erstellt ein neues Verzeichnis mit den gegebenen Pfad.

      \return true wenn erfolgreich
     */
    static bool create_directory(const std::string &path);


private:
    //! Der aktuelle Verzeichnispfad		
    std::string current_directory;
#ifdef __linux__

    //! Vergleicht die Namen von zwei Verzeichniseinträgen
    static bool compare(dirent entry_one, dirent entry_two);

    // Vergleich die Namen von zwei Verzeichniseintraegen anhand der Laenge
    static bool compare_length(dirent entry_one, dirent entry_two);


#else
    //! Vergleicht die Namen von zwei Verzeichniseinträgen
    static bool compare(const WIN32_FIND_DATA &entry_one, const WIN32_FIND_DATA &entry_two);

    //! Vergleicht die Namen von zwei Verzeichniseinträgen anhand der Laenge
    static bool compare_length(const WIN32_FIND_DATA &entry_one, const WIN32_FIND_DATA &entry_two);

    //! Macht ein Verzeichnis Windowskompatible
    /*!
      Window möchte unbedingt einen tollen * haben, so
      kann z.B. C:\ nur ausgelesen werden wenn dem
      OS C:\* übergeben wird. Passt den übergebenen C String an.

      @param[in] path Das Orignalverzeichnis
      @param[out] out Das neue Verzeichnis

     */
    void make_windows_compatible(const TCHAR *path, std::string& out);
#endif


};

#endif // __DIRECTORY_H
