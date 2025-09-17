// Provides access to system resources and compiler functions


use parking_lot::Mutex;

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::JsValue;

use chrono::{
    DateTime, 
    Datelike, 
    Local
};

use std::{
    collections::HashMap,
    path::{Path, PathBuf},
    sync::OnceLock,
    cell::OnceCell
};

use typst::{
    diag::{EcoString, FileError, FileResult},
    foundations::{Bytes, Datetime},
    model::Document,
    syntax::{FileId, Source, VirtualPath},
    text::{Font, FontBook},
    utils::LazyHash,
    Library, World,
};

#[cfg(target_arch = "wasm32")]
type SystemWorldError = JsValue;

#[cfg(not(target_arch = "wasm32"))]
type SystemWorldError = anyhow::Error;


pub struct FileEntry {
    bytes: OnceCell<Bytes>,
    source: Source,
}

impl FileEntry {
    pub fn new(id: FileId, text: String) -> Self {
        Self {
            bytes: OnceCell::new(),
            source: Source::new(id, text),
        }
    }

    pub fn source(&self) -> Source {
        self.source.clone()
    }

    pub fn bytes(&self) -> Bytes {
        self.bytes
            .get_or_init(|| Bytes::from(self.source.text().as_bytes()))
            .clone()
    }
}

pub struct EmbeddedFonts {
    pub regular: &'static [u8],
    pub bold: &'static [u8],
    pub italic: &'static [u8],
    pub bold_italic: &'static [u8],
    pub extra_light: &'static [u8],
    pub extra_light_italic: &'static [u8],
    pub light: &'static [u8],
    pub light_italic: &'static [u8],
    pub medium: &'static [u8],
    pub medium_italic: &'static [u8],
    pub semi_bold: &'static [u8],
    pub semi_bold_italic: &'static [u8],
    pub thin: &'static [u8],
    pub thin_italic: &'static [u8],
}
impl EmbeddedFonts {
    pub fn all_fonts(&self) -> Vec<&'static [u8]> {
        vec![
            self.regular,
            self.bold,
            self.italic,
            self.bold_italic,
            self.extra_light,
            self.extra_light_italic,
            self.light,
            self.light_italic,
            self.medium,
            self.medium_italic,
            self.semi_bold,
            self.semi_bold_italic,
            self.thin,
            self.thin_italic,
        ]
    }
}

pub fn get_embedded_fonts() -> EmbeddedFonts {
    EmbeddedFonts {
        regular: include_bytes!("../fonts/IBMPlexSans-Regular.ttf"),
        bold: include_bytes!("../fonts/IBMPlexSans-Bold.ttf"),
        italic: include_bytes!("../fonts/IBMPlexSans-Italic.ttf"),
        bold_italic: include_bytes!("../fonts/IBMPlexSans-BoldItalic.ttf"),
        extra_light: include_bytes!("../fonts/IBMPlexSans-ExtraLight.ttf"),
        extra_light_italic: include_bytes!("../fonts/IBMPlexSans-ExtraLightItalic.ttf"),
        light: include_bytes!("../fonts/IBMPlexSans-Light.ttf"),
        light_italic: include_bytes!("../fonts/IBMPlexSans-LightItalic.ttf"),
        medium: include_bytes!("../fonts/IBMPlexSans-Medium.ttf"),
        medium_italic: include_bytes!("../fonts/IBMPlexSans-MediumItalic.ttf"),
        semi_bold: include_bytes!("../fonts/IBMPlexSans-SemiBold.ttf"),
        semi_bold_italic: include_bytes!("../fonts/IBMPlexSans-SemiBoldItalic.ttf"),
        thin: include_bytes!("../fonts/IBMPlexSans-Thin.ttf"),
        thin_italic: include_bytes!("../fonts/IBMPlexSans-ThinItalic.ttf"),
    }
}

pub struct SystemWorld {
    /// The root relative to which absolute paths are resolved.
    root: PathBuf,
    /// The input source.
    main: FileId,
    /// Typst's standard library.
    library: LazyHash<Library>,
    /// Metadata about discovered fonts.
    book: LazyHash<FontBook>,
    /// Storage of fonts
    fonts: Vec<Font>,

    files: Mutex<HashMap<FileId, FileEntry>>,
    now: OnceLock<DateTime<Local>>,

}

impl SystemWorld {
    pub fn new(root: String) -> Result<SystemWorld, SystemWorldError> {

        let (book, fonts) = Self::start_embedded_fonts();
        
        Ok(Self {
            root: PathBuf::from(root),
            main: FileId::new(None, VirtualPath::new("")),
            library: LazyHash::new(Library::default()),
            book: LazyHash::new(book),
            fonts,
            files: Mutex::default(),
            now: OnceLock::new(),
        })
    }

    pub fn compile(&mut self, text: String, vpath: String, report_logo: Option<Vec<u8>>) -> Result<Document, SystemWorldError> {

        // Check if this adds too much overhead
        if let Some(data) = report_logo {
            self.add_logo(data)?;
        } else {
            self.add_logo(SystemWorld::default_logo())?;
        }

        self.main = FileId::new(None, VirtualPath::new(vpath));

        self.files
            .get_mut()
            .insert(self.main, FileEntry::new(self.main, text));

        let document = typst::compile(self)
            .output
            .map_err(|_| {

                #[cfg(target_arch = "wasm32")]
                let err = JsValue::from("failed to compile document");
                #[cfg(not(target_arch = "wasm32"))]
                let err = anyhow::anyhow!("failed to compile document");

                // Deallocate sensitive data from virtual path memory on compilation error
                self.reset();

                err
        });

        // Deallocate sensitive data from virtual path memory after compilation of the document
        self.reset();

        document
    }

    pub fn default_logo() -> Vec<u8> {
        include_bytes!("../logos/vidrl.png").to_vec()
    }
    pub fn certificate_logo() -> Vec<u8> {
        include_bytes!("../logos/cert.png").to_vec()
    }
    
    pub fn add_logo(&mut self, data: Vec<u8>) -> Result<(), SystemWorldError> {

        let logo_id = FileId::new(
            None, 
            VirtualPath::new("logo.png")
        );

        let logo_bytes = Bytes::from(data);

        let logo_entry = FileEntry {
            bytes: OnceCell::new(),
            source: Source::new(logo_id, String::new()), // No text source needed for binary data
        };

        if let Err(_) = logo_entry.bytes.set(logo_bytes) {
            #[cfg(target_arch = "wasm32")]
                let err = JsValue::from("failed to set logo bytes in OnceCell");
                #[cfg(not(target_arch = "wasm32"))]
                let err = anyhow::anyhow!("failed to set logo bytes in OnceCell");
                return Err(err )
        };

        self.files.lock().insert(logo_id, logo_entry);

        Ok(())
    }

    fn reset(&mut self) {
        self.files = Mutex::new(HashMap::new()); // fully deallocate files from memory
        self.main = FileId::new(None, VirtualPath::new("")); // reset and deallocate main file
        self.now.take();
    }

    fn read_file(&self, path: &Path) -> FileResult<String> {
        std::fs::read_to_string(path).map_err(|e| match e.kind() {
            std::io::ErrorKind::NotFound => FileError::NotFound(path.to_path_buf()),
            std::io::ErrorKind::PermissionDenied => FileError::AccessDenied,
            _ => FileError::Other(Some(EcoString::from(format!("{}", e)))),
        })
    }


    fn file_entry<F, T>(&self, id: FileId, f: F) -> FileResult<T>
    where
        F: FnOnce(&mut FileEntry) -> T,
    {
        let mut map = self.files.lock();
        if !map.contains_key(&id) {
            let path = self.root.clone();    
            let resolved_path = id
                .vpath()
                .resolve(&path)
                .ok_or(FileError::AccessDenied)?;

            let text = self.read_file(&resolved_path)?;

            map.insert(id, FileEntry::new(id, text));
        }
        Ok(f(map.get_mut(&id).unwrap()))
    }

    fn start_embedded_fonts() -> (FontBook, Vec<Font>) {
        let mut book = FontBook::new();
        let mut fonts = Vec::new();

        for data in typst_assets::fonts() {
            let buffer = Bytes::from_static(data);
            for font in Font::iter(buffer) {
                book.push(font.info().clone());
                fonts.push(font);
            }
        }
        let embedded_fonts = get_embedded_fonts();
        
        for data in embedded_fonts.all_fonts() {
            let buffer = Bytes::from_static(data);
            for font in Font::iter(buffer) {
                book.push(font.info().clone());
                fonts.push(font);
            }
        }

        return (book, fonts);
    }
}

impl World for SystemWorld {
    fn library(&self) -> &LazyHash<Library> {
        &self.library
    }
    fn book(&self) -> &LazyHash<FontBook> {
        &self.book
    }
    fn main(&self) -> FileId {
        self.main
    }
    fn source(&self, id: FileId) -> FileResult<Source> {
        self.file_entry(id, |f| f.source())
    }
    fn file(&self, id: FileId) -> FileResult<Bytes> {
        self.file_entry(id, |f| f.bytes())
    }
    fn font(&self, index: usize) -> Option<Font> {
        Some(self.fonts[index].clone())
    }
    fn today(&self, offset: Option<i64>) -> Option<Datetime> {
        let now = self.now.get_or_init(chrono::Local::now);

        let naive = match offset {
            None => now.naive_local(),
            Some(o) => now.naive_utc() + chrono::Duration::hours(o),
        };

        Datetime::from_ymd(
            naive.year(),
            naive.month().try_into().ok()?,
            naive.day().try_into().ok()?,
        )
    }
}